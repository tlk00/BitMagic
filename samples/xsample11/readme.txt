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
for volume, and sparse_vector_float's for the exchange rates and percent changes,
Optimizes the sparse_vectors, and freezes them as no new information will be added.

The program then gathers information on how much memory the sparse_vector's
use via calc_stat, and then serializes the vectors to check their serialized size on
disk.

The use case illustrated shows the pct_changes being searched with a scanner for
whenever the EUR to USD percent change is above 1%, and when the USD to JPY 
percent change is below -1%, meaning that the US dollar likely fluctuated and decreased
in value whenever the two match.

Some other possible use cases for this financial data:
- Searching for when the USD increased in value by 1%
- Searching for when 1 USD was worth between 90 and 100 JPY
- Searching for if 1 EUR was worth less than 1 USD in the timespan

Data Shown in sample:
Size of data structures in memory in bytes

	    std::vector<float> | 399876
     std::vector<unsigned int> | 399876
      std::vector<std::string> | 2399256
--------------------------------------------------------
                eur_day(dates) | 284350
--------------------------------------------------------
                      eur_open | 375760
                      eur_high | 375236
                       eur_low | 375760
                eur_pct_change | 536960
                     eur_close | 375236
--------------------------------------------------------
                      jpy_open | 370448
                      jpy_high | 369380
                       jpy_low | 370448
                jpy_pct_change | 536576
                     jpy_close | 369380

Serialized size of the data in the sparse_vectors

Dates Serialized Size:          82551 bytes
--------------------------------------------------------
EUR Open Serialized Size:       229206 bytes
EUR High Serialized Size:       228614 bytes
EUR Low Serialized Size:        229206 bytes
EUR Pct Change Serialized Size: 373193 bytes
EUR Close Serialized Size:      228614 bytes
EUR Volume Serialized Size:     194193 bytes
--------------------------------------------------------
JPY Open Serialized Size:       233855 bytes
JPY High Serialized Size:       232499 bytes
JPY Low Serialized Size:        233855 bytes
JPY Pct Change Serialized Size: 371879 bytes
JPY Close Serialized Size:      232499 bytes
JPY Volume Serialized Size:     190561 bytes
--------------------------------------------------------
COMBINED TOTAL SERIALIZED SIZE: 3060725 bytes




