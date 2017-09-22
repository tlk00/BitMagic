#ifndef BMTIMER__H__INCLUDED__
#define BMTIMER__H__INCLUDED__
/*
Copyright(c) 2002-2017 Anatoliy Kuznetsov(anatoliy_kuznetsov at yahoo.com)

Permission is hereby granted, free of charge, to any person 
obtaining a copy of this software and associated documentation 
files (the "Software"), to deal in the Software without restriction, 
including without limitation the rights to use, copy, modify, merge, 
publish, distribute, sublicense, and/or sell copies of the Software, 
and to permit persons to whom the Software is furnished to do so, 
subject to the following conditions:

The above copyright notice and this permission notice shall be included 
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, 
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, 
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR 
OTHER DEALINGS IN THE SOFTWARE.

For more information please visit:  http://bmagic.sourceforge.net

*/

// BitMagic timing functions (internal header)

#include <iostream>
#include <map>
#include <chrono>

namespace bm
{


/// Utility class to collect performance measurements and statistics.
///
class chrono_taker
{
public:
    /// collected statistics
    ///
    struct statistics
    {
        std::chrono::duration<double, std::milli>  duration;
        unsigned                                   repeats;
        
        statistics() : repeats(1) {}
        
        statistics(std::chrono::duration<double, std::milli> d, unsigned r)
        : duration(d), repeats(r)
        {}
    };
    
    /// test name to duration map
    ///
    typedef std::map<std::string, statistics > duration_map_type;

public:
    chrono_taker(const std::string name,
                unsigned repeats = 1,
                duration_map_type* dmap = 0)
    : name_(name),
      repeats_(repeats),
      dmap_(dmap),
      is_stopped_(false)
    {
        start_ = std::chrono::steady_clock::now();
    }
    
    ~chrono_taker()
    {
        try
        {
            if (!is_stopped_)
            {
                stop();
            }
        }
        catch(...)
        {}
    }
    
    
    void stop(bool silent=false)
    {
        finish_ = std::chrono::steady_clock::now();
        auto diff = finish_ - start_;
        if (dmap_)
        {
            (*dmap_)[name_] = statistics(diff, repeats_);
        }
        else // report the measurements
        {
            if (!silent)
                std::cout << name_ << ": "
                          << std::chrono::duration <double, std::milli> (diff).count()
                          << std::endl;
        }
        is_stopped_ = true;
    }
    
    static void print_duration_map(const duration_map_type& dmap)
    {
        duration_map_type::const_iterator it = dmap.begin();
        duration_map_type::const_iterator it_end = dmap.end();
        
        for (;it != it_end; ++it)
        {
            std::cout << it->first << "; " << it->second.duration.count() << std::endl;
        }
    }
    
    

    chrono_taker(const chrono_taker&) = delete;
    chrono_taker & operator=(const chrono_taker) = delete;
    
protected:
    std::string                                        name_;
    std::chrono::time_point<std::chrono::steady_clock> start_;
    std::chrono::time_point<std::chrono::steady_clock> finish_;
    unsigned                                           repeats_;
    duration_map_type*                                 dmap_;
    bool                                               is_stopped_;
};


} // namespace

#endif
