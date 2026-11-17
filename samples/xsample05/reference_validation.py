"""Fail-fast xsample05 reference validation with retained, external case logs."""
import concurrent.futures
import json
import os
import signal
import subprocess
import threading
import time
from pathlib import Path


def validate(args, inputs, selected, check_case):
    """Validate selected cases; terminate active readers on the first failure or signal."""
    root = args.out.resolve()
    if root == args.reference_dir or args.reference_dir in root.parents:
        raise ValueError('--out must be outside the immutable reference directory')
    root.mkdir(parents=True, exist_ok=False)
    stop = threading.Event()
    lock = threading.RLock()
    children = set()
    first = []
    results = []

    def abort(message):
        """Preserve the first failure and stop every registered reader."""
        with lock:
            if not first:
                first.append(message)
            stop.set()
            active = list(children)
        for child in active:
            if child.poll() is None:
                try:
                    os.killpg(child.pid, signal.SIGTERM)
                except ProcessLookupError:
                    pass
        deadline = time.monotonic()+2
        for child in active:
            try:
                child.wait(timeout=max(.01, deadline-time.monotonic()))
            except subprocess.TimeoutExpired:
                try:
                    os.killpg(child.pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
                child.wait()
            finally:
                try:
                    os.killpg(child.pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass

    def interrupted(signum, frame):
        """Convert external interruption to the same cancellation and reporting path."""
        abort('Interrupted by ' + signal.Signals(signum).name)

    previous = {sig: signal.signal(sig, interrupted) for sig in (signal.SIGINT, signal.SIGTERM, signal.SIGHUP)}

    def task(case):
        """Keep command, output and timing for one case, including preflight failures."""
        name = case[0]
        if stop.is_set():
            return dict(case=name, status='cancelled', seconds=0)
        work = root / name
        work.mkdir()
        log = work / 'validation.log'
        started = time.monotonic()
        try:
            check_case(case, args.reference_dir, inputs, stop.is_set)
            source = inputs['remap' if case[2] else 'plain']['path']
            command = [str(args.bin), '-svin', str(args.reference_dir / 'cases' / name / 'reference.sv'), '-verify-sv', source]
            (work / 'command.json').write_text(json.dumps(command, indent=2)+'\n')
            with log.open('w') as output:
                with lock:
                    if stop.is_set():
                        return dict(case=name, status='cancelled', seconds=0)
                    child = subprocess.Popen(command, stdout=output, stderr=subprocess.STDOUT, cwd=work, start_new_session=True,
                                             env=dict(os.environ, TMPDIR=str(work), TMP=str(work), TEMP=str(work)))
                    children.add(child)
                try:
                    code = child.wait()
                finally:
                    try:
                        os.killpg(child.pid, signal.SIGKILL)
                    except ProcessLookupError:
                        pass
                    with lock:
                        children.discard(child)
            if code:
                reason = signal.Signals(-code).name if code < 0 else 'exit ' + str(code)
                raise RuntimeError(f'{reason} | command {command} | log {log}')
            return dict(case=name, status='pass', seconds=time.monotonic()-started, log=str(log))
        except Exception as error:
            if stop.is_set():
                return dict(case=name, status='cancelled', seconds=time.monotonic()-started, log=str(log))
            message = f'{name}: {error} | log {log}'
            with log.open('a') as output:
                output.write('\nVALIDATION ERROR: '+message+'\n')
            abort(message)
            return dict(case=name, status='failed', seconds=time.monotonic()-started, log=str(log), error=str(error))

    try:
        with concurrent.futures.ThreadPoolExecutor(max_workers=args.j) as pool:
            futures = [pool.submit(task, case) for case in selected]
            for future in concurrent.futures.as_completed(futures):
                try:
                    result = future.result()
                except Exception as error:
                    abort(str(error))
                    result = dict(case='runner', status='failed', error=str(error))
                results.append(result)
                passed = sum(r['status'] == 'pass' for r in results)
                report = dict(total=len(selected), passed=passed, results=results, failure=first[0] if first else None)
                temporary = root / 'progress.json.tmp'
                temporary.write_text(json.dumps(report, indent=2)+'\n')
                temporary.replace(root / 'progress.json')
                print(f"{result['status'].upper()} {result['case']} | {passed}/{len(selected)} passed | {result.get('log', root)}", flush=True)
    finally:
        for sig, handler in previous.items():
            signal.signal(sig, handler)
    if first:
        raise RuntimeError(first[0])
    print(f'validate complete: {len(selected)} cases | results {root}', flush=True)
