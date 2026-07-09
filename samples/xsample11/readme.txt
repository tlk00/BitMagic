Example to illustrate dataset scrolling using different de-serialization
techniques with succinct vectors.

Example to illustrate how to use a sparse_vector_float in a real-world
setting

https://bitmagic.io/bm-9.1.0


How to build and run:

Unzip EURUSD_H1.csv.zip
Unzip USDJPY_H1.csv.zip

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

The program then gathers information on how much memory the sparse_vector_float's
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


