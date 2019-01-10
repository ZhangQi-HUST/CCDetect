# CCDetect
A new chessboard corner detection algorithm with simple thresholding

This project is the implementation of our chessboard corner detection algorithm.
We used MATLAB to do most of the work, such as image reading, filtering, and saving. Due to that MATLAB is quite slow to process images pixel by pixel, we adopt the MATLAB/C++ hybrid programming. The new corner response function is written in C++. We only implemented a circular sampling with 16 samples.

There are three cpp files uploaded. 
- mexLCFP.cpp  input parameter (image, sigma)
- mexLCFPn.cpp  input parameter (image, sigma, radius)
- mexLCFPnt.cpp  input parameter (image, threshold, radius)

The output of these function are (response, mask).
response is the response map.
mask store the local maxima. You can get all position of local maximum by 
``` 
[y,x]=find(mask>=4)
```

These cpp files need to be compiled with `mex` command.

Please reference our work
Zhang, Q. and C. Xiong. A New Chessboard Corner Detection Algorithm with Simple Thresholding. 2017. Cham: Springer International Publishing.
