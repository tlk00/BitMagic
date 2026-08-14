# rscsample05: Collaborative XOR compression

This example serializes a data-frame-like group of RSC and regular sparse
vectors with collaborative XOR compression. Similar bit-planes from other
columns act as references, reducing the combined serialized size.

The program compares ordinary independent serialization with the
collaborative form. It builds a shared reference list in reverse column order,
computes one XOR similarity model, attaches it to serializers of different
vector types and stores all vector BLOBs behind a simple size header.

Collaborative deserialization is a two-pass operation: first restore every
vector's structure and build the reverse-order reference list, then decode the
data without clearing those structures. Serializers and deserializers are
detached from their non-owning reference lists before the lists leave scope.
