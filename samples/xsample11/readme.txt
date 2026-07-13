Example to illustrate dataset scrolling using different de-serialization
techniques with succinct vectors.

Example to illustrate how to use a sparse_vector_float in a real-world
setting

https://bitmagic.io/bm-9.1.0


How to build and run:

Unzip EURUSD_H1.csv.zip
Unzip USDJPY_H1.csv.zip

Make sure these files are in the same directory as where you are running the xsample11 from

make

How to run:
-------------
./xsample11


Application notes:
------------------

Example reads the data in the EURUSD_H1.csv and the USDJPY_H1.csv 
The data in the .csv's is the exchange rate information between EUR to USD and 
USD to JPY hourly from June 24, 2010 5:00am to July 7, 2026 3:00pm, (15:00)
This information includes each hours open, close, high, low rates, and volume.

This program then calculates each hour's percent change based on the close, and
stores all of the information in a str_sparse_vector for times, sparse_vector
for volume, and sparse_vector_float's for the exchange rates, and a sparse_vector of ints for,
pct changes. 
Optimizes the sparse_vectors, and freezes them as no new information will be added.

The program then gathers information on how much memory the sparse_vector's
use via calc_stat, and then serializes the vectors to check their serialized size on
disk.

One use case illustrated shows the pct_changes being searched with a scanner for
whenever the EUR to USD percent change is above 1%, and when the USD to JPY 
percent change is below -1%, meaning that the US dollar likely fluctuated and decreased
in value whenever the two match.

Another use case illustrated shows the number of times the data appeared in a range near
certain values using a scanner, When the USD to JPY rate was 1 USD was near 100 JPY and 103 JPY
to see if being near 100 JPY Has some psychological effect that made the rate stay 1:100
for longer than a random other rate which is not as significant.

Some other possible use cases for this financial data:
- Searching for when the USD increased in value by 1%
- Searching for when 1 USD was worth between 90 and 100 JPY
- Searching for if 1 EUR was worth less than 1 USD in the timespan

Data Shown in sample:
Size of data structures in memory in bytes

	    std::vector<float> | 399876
     std::vector<unsigned int> | 399876
      std::vector<std::string> | 1899411 (Not including std::string overhead)
--------------------------------------------------------
                eur_day(dates) | 284350
       jpy_day(dates remapped) | 238762
--------------------------------------------------------
                      eur_open | 375760
                      eur_high | 375236
                       eur_low | 375760
                eur_pct_change | 191864
                     eur_close | 375236
--------------------------------------------------------
                      jpy_open | 370448
                      jpy_high | 369380
                       jpy_low | 370448
                jpy_pct_change | 194416
                     jpy_close | 369380
--------------------------------------------------------
Total EUR sparse_vector memory usage: 1992860 bytes
Total JPY sparse_vector memory usage: 1971288 bytes

Total memory usage: 4248498 bytes
Total memory usage using std::vector's: 6697923 bytes


Serialized size of the data in the sparse_vectors

Dates Serialized Size:          82551 bytes
Remapped Dates Serialized Size: 84990 bytes
--------------------------------------------------------
EUR Open Serialized Size:       229206 bytes
EUR High Serialized Size:       228614 bytes
EUR Low Serialized Size:        229206 bytes
EUR Pct Change Serialized Size: 114418 bytes
EUR Close Serialized Size:      228614 bytes
EUR Volume Serialized Size:     194193 bytes
--------------------------------------------------------
JPY Open Serialized Size:       233855 bytes
JPY High Serialized Size:       232499 bytes
JPY Low Serialized Size:        233855 bytes
JPY Pct Change Serialized Size: 115502 bytes
JPY Close Serialized Size:      232499 bytes
JPY Volume Serialized Size:     190561 bytes
--------------------------------------------------------
COMBINED TOTAL SERIALIZED SIZE: 2545573 bytes


Time to run a single float range search with scanner: 0.108709 to 0.08125 ms

