#ifndef BM_CMD_ARGS_H__INCLUDED__
#define BM_CMD_ARGS_H__INCLUDED__

inline
void show_help()
{
    std::cerr
        << "BitMagic DNA k-mer build and count (c) 2020" << std::endl
        << "-fa   file-name            -- input FASTA file" << std::endl
        << "-k    size                 -- k-mer size (4,8,16,..24) " << std::endl
        << "-kd   file-name            -- k-mer dictionary file (output)"  << std::endl
        << "-kdf  file-name            -- k-mer dictionary for frequent k-mers (see -fpc)" << std::endl
        << "-kdr  file-name            -- k-mer dictionary clean from frequent(over-represented) k-mers" << std::endl
        << "-fpc  <int>                -- frequent k-mers to exclude (see -kdr) percents " << std::endl
        << "-kdc  file-name            -- k-mer counts file (output)" << std::endl
        << "-j    number-of-threads    -- number of parallel jobs to run" << std::endl
        << "-kh   file-name            -- k-mer counts distibution histogram (TSV)" << std::endl
        << "-diag                      -- run diagnostics"  << std::endl
        << "-timing                    -- collect timings"  << std::endl
      ;
}

// cmd-line arguments parser
//
inline
int parse_args(int argc, char *argv[])
{
    for (int i = 1; i < argc; ++i)
    {
        std::string arg = argv[i];
        if ((arg == "-h") || (arg == "--help"))
        {
            show_help();
            return 0;
        }

        if (arg == "-fa" || arg == "--fa")
        {
            if (i + 1 < argc)
            {
                ifa_name = argv[++i];
            }
            else
            {
                std::cerr << "Error: -fa requires file name" << std::endl;
                return 1;
            }
            continue;
        }
        if (arg == "-kd" || arg == "--kd")
        {
            if (i + 1 < argc)
            {
                ikd_name = argv[++i];
            }
            else
            {
                std::cerr << "Error: -kd requires file name" << std::endl;
                return 1;
            }
            continue;
        }
        if (arg == "-kdc" || arg == "--kdc")
        {
            if (i + 1 < argc)
            {
                ikd_counts_name = argv[++i];
            }
            else
            {
                std::cerr << "Error: -kdc requires file name" << std::endl;
                return 1;
            }
            continue;
        }
        if (arg == "-kdr" || arg == "--kdr")
        {
            if (i + 1 < argc)
            {
                ikd_rep_name = argv[++i];
            }
            else
            {
                std::cerr << "Error: -kdr requires file name" << std::endl;
                return 1;
            }
            continue;
        }
        if (arg == "-kdf" || arg == "--kdf")
        {
            if (i + 1 < argc)
            {
                ikd_freq_name = argv[++i];
            }
            else
            {
                std::cerr << "Error: -kdf requires file name" << std::endl;
                return 1;
            }
            continue;
        }

        if (arg == "-kh" || arg == "--kh")
        {
            if (i + 1 < argc)
            {
                kh_name = argv[++i];
            }
            else
            {
                std::cerr << "Error: -kh requires file name" << std::endl;
                return 1;
            }
            continue;
        }

        if (arg == "-k" || arg == "--k")
        {
            if (i + 1 < argc)
            {
                std::string ksize_str = argv[++i];
                if (ksize_str.empty())
                {
                    std::cerr << "Error: -k requires size" << std::endl;
                    return 1;
                }
                ik_size = (unsigned)std::stoi(ksize_str);
                if (ik_size < 2)
                {
                    ik_size = 8;
                }
                if (ik_size > 24)
                {
                    std::cerr << "Error: unsupported k-mer size (too big)" << std::endl;
                    return 1;
                }
            }
            else
            {
                std::cerr << "Error: -k requires size" << std::endl;
                return 1;
            }
            continue;
        }

        if (arg == "-fpc" || arg == "--fpc")
        {
            if (i + 1 < argc)
            {
                f_percent = unsigned(::atoi(argv[++i]));
                if (f_percent == 0 || f_percent >= 100)
                {
                    std::cerr << "WARNING: -fpc value incorrect " << f_percent;
                    f_percent = 5;
                    std::cerr << " assumed as " << f_percent << endl;
                }
            }
            else
            {
                std::cerr << "Error: -fpc requires integer " << std::endl;
                return 1;
            }
            continue;
        }

        if (arg == "-j" || arg == "--j")
        {
            if (i + 1 < argc)
            {
                parallel_jobs = unsigned(::atoi(argv[++i]));
            }
            else
            {
                std::cerr << "Error: -j requires number of jobs" << std::endl;
                return 1;
            }
            continue;
        }


        if (arg == "-diag" || arg == "--diag" || arg == "-d" || arg == "--d")
            is_diag = true;
        if (arg == "-timing" || arg == "--timing" || arg == "-t" || arg == "--t")
            is_timing = true;
        if (arg == "-bench" || arg == "--bench" || arg == "-b" || arg == "--b")
            is_bench = true;


    } // for i
    return 0;
}


#endif

