isONclust2
==========

A minimal example
------------------

```bash
# sort reads and write out batches:
isONclust2 sort -B 50000 -v ens500.fq
# initial clustering of individual batches:
isONclust2 cluster -v -l isONclust2_batches/batches/isONbatch_0.cer -o b0.cer
isONclust2 cluster -v -l isONclust2_batches/batches/isONbatch_1.cer -o b1.cer
isONclust2 cluster -v -l isONclust2_batches/batches/isONbatch_2.cer -o b1.cer
# merge cluster batches:
isONclust2 cluster -v -l b0.cer -r b1.cer -o b_0_1.cer
isONclust2 cluster -v -l b_0_1.cer -r b2.cer -o b_0_1_2.cer
# dump final results:
isONclust2 dump -v -o results b_0_1_2.cer
```

Building from source
--------------------

- Clone the repository: `git clone --recursive https://github.com/nanoporetech/isONclust2.git`
- Make sure you have cmake v3.1 or later.
- Issue `cd isONclust2; mkdir build; cd build; cmake ..; make -j`
