# svsample07a

Example for bit-transposed succinct vector.

Example illustrates use of a scanner to do a quick search in the unordered bit-transposed vector.
Results bit-vector is inspected using bm::interval_enumerator<> to find patterns:

1. unique element
2. elements co-located in the search vector together 
3. elements dispersed in the search vector (multiple occasions).

## Sample run

./svsample07a 

> Sparse vector:
> 
> 0,4294967295,17,17,5,18,178,178,17,0,4294967295,
> 
> Value = 0 is not colocated
> 
> Value = 4294967295 is not colocated
> 
> Value = 17 is not colocated
> 
> Value = 5 is unique
> 
> Value = 18 is unique
> 
> Value = 178 is colocated
