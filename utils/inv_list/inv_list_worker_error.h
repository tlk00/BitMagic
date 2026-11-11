#ifndef INV_LIST_WORKER_ERROR_H_INCLUDED
#define INV_LIST_WORKER_ERROR_H_INCLUDED
/*
Copyright(c) 2026 Anatoliy Kuznetsov
Licensed under the Apache License, Version 2.0.
https://www.apache.org/licenses/LICENSE-2.0
*/
/** @file inv_list_worker_error.h
    @brief Main-thread reporting for completed inv_list worker failures.
    @note Internal application helper; call only after acquiring task completion.
*/
#include <exception>
#include <stdexcept>
#include <string>

/** Report a completed worker failure with operation context and original cause.
    @param error Captured worker exception; may be empty for pool-level failures.
    @param code Completed task error code; zero indicates success.
    @param context Main-thread phase, file and vector information.
    @pre The task's completion flag was observed with acquire semantics.
    @throws std::runtime_error If either failure signal is set. Standard exception
    messages are retained; non-standard exceptions and bare error codes are named.
*/
inline
void inv_list_check_worker(std::exception_ptr error, int code,
                           const std::string& context)
{
    if (error)
    {
        try
        {
            std::rethrow_exception(error);
        }
        catch (const std::exception& cause)
        {
            throw std::runtime_error(context + ": " + cause.what());
        }
        catch (...)
        {
            throw std::runtime_error(context + ": non-standard worker exception");
        }
    }
    if (code)
    {
        throw std::runtime_error(context + ": task error code=" + std::to_string(code));
    }
}
#endif
