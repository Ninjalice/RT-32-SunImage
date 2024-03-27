# RT-32
---

## TODO 

- [x] Fix UTC time data tipe from jd and decimal hour to decimal seconds.
- [x] Fix interpolaion of antena movement for merging with real data.
- [x] Create basic image with data points.
- [ ] Apply convolution to image.
  - [x] Apply basic convolution kernel. (Bad result)
  - [x] Apply manual non zero kernel
  - [ ] Try diferent kernels
- [ ]  Apply scipy.signal savgol_filter
- [ ] Remove calibration points and calibrate it correctly
- [ ] Try to apply gaussian blurr and noise reduction
- [ ] Use real antena movements rather than theoretical positions.
- [ ] Generalise code for all bands and polarizations
  

## FUTURE 

- [ ] Think about repository and project hosting.
- [ ] Create React control webpage.
- [ ] Think of automatic band change in Arduino.
- [ ] Think of automatic uploading of antena positions for plc.
