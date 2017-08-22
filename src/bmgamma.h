/*
Copyright(c) 2002-2009 Anatoliy Kuznetsov(anatoliy_kuznetsov at yahoo.com)

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

#ifndef BMGAMMAENC__H__INCLUDED__
#define BMGAMMAENC__H__INCLUDED__


namespace bm
{


/**
    Elias Gamma decoder
*/
template<typename T, typename TBitIO>
class gamma_decoder
{
public:
    gamma_decoder(TBitIO& bin) : bin_(bin) 
    {}
    
    /**
        Start encoding sequence
    */
    void start()
    {}
    
    /**
        Stop decoding sequence
    */
    void stop()
    {}
    
    /**
        Decode word
    */
    T operator()(void)
    {
        unsigned l = bin_.eat_zero_bits();
        bin_.get_bit(); // get border bit
        T current = 0;
        for (unsigned i = 0; i < l; ++i)
        {
            if (bin_.get_bit())
            {
                current += 1 << i;
            }
        }
        current |= (1 << l);
        return current;
    }
private:
    gamma_decoder(const gamma_decoder&);
    gamma_decoder& operator=(const gamma_decoder&);
private:
    TBitIO&  bin_;
};



} // bm

#endif

