# PCA-PMI
Codes for "Part mutual information for quantifying direct associations in networks"
https://doi.org/10.1073/pnas.1522586113

"PCA" is nothing to do with Principle Component Analysis. Here, it is PC algorithm ([Peter Spirtes, Clark Glymour](https://link.springer.com/book/10.1007/978-1-4612-2748-9)), a network structure inference algorithm. PMI is part mutual information, a new criteria to estimate condition independence.

## How to use the codes

All of the methods we implemented and the methods we want to compare with are in the folder `/lib`.

Please  run `make.m` to add the path into Matlab Search path temporarily.

If you want to add the path permanently, please run `make -p` . You can remove the path from "Set Path" ![path](SetPath.png)bottom.

- `cmi.m`  
  - Computing Conditional Mutual Information with bin method
- `pmi.m`
  - Computing Part Mutual Information with bin method
- `pmiguass.m`
  - Computing Part Mutual Information with Gaussian Distribution Approximation 
- `dcorr.m`
  - Computing Distance Correlation
- `pdcor`
  - Computing Partial Distance Correlation
- `graphicalLasso.m`
  - Network Structure Inference with graphicalLasso
- `pca_pmi.m`
  - Network Structure Inference with PC algorithm combining with Part Mutual Information
- `kpca_pmi.m`
  - Network Structure Inference with PC algorithm combining with Kernelized Part Mutual Information
- `pca_cmi.m`
  - Network Structure Inference with PC algorithm combining with Conditional Mutual Information
- `pcapcc.m`
  - Network Structure Inference with PC algorithm combining with Pearson Correlation 

## Example

Please run `configpath.m` in the folder `example` first.

Please read `example/README.md` to get all the descriptions of codes in `example` .

## Running Environment

* Windows, Unix/Linux, Mac OS
* Matlab (>2013b)

## Citation
@article{zhao2016part,  
  title={Part mutual information for quantifying direct associations in networks},  
  author={Zhao, Juan and Zhou, Yiwei and Zhang, Xiujun and Chen, Luonan},  
  journal={Proceedings of the National Academy of Sciences},
  volume={113},  
  number={18},  
  pages={5130--5135},  
  year={2016},  
  publisher={National Acad Sciences}  
}

## License

Apache License 2.0