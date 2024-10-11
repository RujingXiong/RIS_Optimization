# Optimal Discrete Beamforming of RIS-Aided Wireless Communications: An Inner Product Maximization Approach
## Description
1. The main branch is for the discrete optimizaion problem with quadratic forms ($w^HRw$) or 2-norms $||Aw||_2$, where the matrix R ( or A) is with rank-1. The approach is called DaS method !!!
2. The For-lowrank branch code is used to achieve the global optimum in maximizing quadratic forms ($w^HRw$) or 2-norms $||Aw||_2$ where the matrix R ( or A) is M, M is not 1. w is the discrete phase configuration vector.
## Reference 1
The main code is associated with the published paper entitled "Optimal Discrete Beamforming of RIS-Aided Wireless Communications: An Inner Product Maximization Approach". If you are interested in our work, please reference it using
```
@inproceedings{xiong2024optimal,
  title={Optimal discrete beamforming of RIS-aided wireless communications: An inner product maximization approach},
  author={Xiong, Rujing and Dong, Xuehui and Mi, Tiebin and Wan, Kai and Qiu, Robert Caiming},
  booktitle={2024 IEEE Wireless Communications and Networking Conference (WCNC)},
  pages={1--6},
  year={2024},
  organization={IEEE}
}
```

## Description 2 
The special code "Alternating_main" in main branch is the Alternating inner product maximization approach, mainly for solving p-norm $||Aw||_p$ (p is 2 within this code) maximization problems, while A is a matrix with rank-M. This approach features high efficiency and is based on the optimal discrete optimization with rank-1 scenarios (DaS method), as stated as above reference 1.

## Reference 2
This code is associated with the paper--
Xiong R, Mi T, Lu J, et al. Optimal Beamforming of RIS-Aided Wireless Communications: An Alternating Inner Product Maximization Approach[J]. arXiv preprint arXiv:2405.06442, 2024. If you are interested in our work, please reference it using
```
@misc{xiong2024optimal,
  title={Optimal Beamforming of RIS-Aided Wireless Communications: An Alternating Inner Product Maximization Approach}, 
      author={Rujing Xiong and Tiebin Mi and Jialong Lu and Ke Yin and Kai Wan and Fuhai Wang and Robert Caiming Qiu},
      year={2024},
      eprint={2405.06442},
      archivePrefix={arXiv},
      primaryClass={cs.IT},
      url={https://arxiv.org/abs/2405.06442}, 
}
```

