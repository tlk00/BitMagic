This example illustrates creation of layout model using succinct data structures.

As a use case example illustrates hypothetical DNA sequences with features 
annotated in sequence coordinates [from...to]. Example illustrates use of
bit-vector and compact (rank-select compressed) bit-transposed sparse vector
to keep feature coordinates and attributes (example uses strand direction to 
illustrate the case). 

Models are visualization friendly, this example uses ASCII art to generate 
pseudo graphics like:

-------------------------------------------------------------------------
ATGTTAGCCCGCGCATATTATATATGTAGCGTATTAAGCGDGGAGATTACCCTTGCATTAGGTTANNNNNNNN
-------------------------------------------------------------------------
<    >....>    <....<
    <.................................................................<
                    >.........>
                ?....?

Legend:
>....> - feature on positive strand
<....< - feature on negative strand
?....? - feature on unknown strand


Another use case this example illustrates: data model slicing.

Slicing targets two main features:
1. Create a range slice of data in genomics coordinates
2. Do "Feature space slicing" - drop vectors represented features we do not need
and save on traffic and compute. Approach typical for all columnar databases.

Slicing example creates a slice of coordinates [from..to] and drops or includes 
features (in this case just strand feature). 

Slicing works with attention fo real span, so if features cross the slicing 
boundaries - slicing boundaries get extended accordingly, so features are 
never cut partial. 

For example slice [5..10] would include all features overlapping this
coordinate boundaries.

-------------------------------------------------------------------------
ATGTTAGCCCGCGCATATTATATATGTAGCGTATTAAGCGDGGAGATTACCCTTGCATTAGGTTANNNNNNNN
-------------------------------------------------------------------------
     >....>
    <.................................................................<


See also:
http://bitmagic.io/gen-layout.html

