# flutter-echo-modeling
This repository presents a Matlab demo for the synthesis of flutter echo. 
## Overview
In this work, we introduce a method to recreate flutter echo as an audio effect. 
The proposed algorithm is based on a feedback structure utilizing velvet noise that aims to synthesize the fluttery components of a reference room impulse response presenting flutter echo. Among these, the repetition time defines the length of the delay line in a feedback filter. The specific spectral properties of the flutter are obtained with a bandpass attenuation filter and a ripple filter, which enhances the harmonic behavior of the sound. Additional temporal shaping of a velvet-noise filter, which processes the output of the feedback loop, is performed based on the properties of the reference flutter
<p align="center">
  <img width="300" src=".\pictures\structure.PNG">
</p>

## Getting started
The project repo contains the demo code `demo.mlx` written in MATLAB Live Editor notebook. Functions and data that are necessary to execute the code are already in the repository. The algorithm reads the reference RIRs stored in `.\audio\RIRs`, extracts the flutter echo parameters, and synthesize the fluttery components based on the extracted parameters. The output of the convolution between the input signal with the synthesized RIR will then be saved in  `.\outputs`.

Required toolboxes:
- signal_toolbox
- wavelet_toolbox

## Author 
**Gloria Dal Santo**   
M.Sc. in Electrical and Electronic Engineering at EPFL  
B.Sc. in Electronic and Communications Engineering at Politecnico di Torino  
[Linkedin](https://www.linkedin.com/in/gloriadalsanto/)  

For any problem, issue or comment please send an email to gloria.dalsanto@outlook.it  
## References
If you would like to use this code, please cite the related research work using the following reference:
```
G. Dal Santo, K. Prawda, and V. Välimäki,, "Flutter Echo Modeling", 
in Proc. 25rd Int. Conf. Digital Audio Effects (DAFx20in22), Vienna, Austria, Sept. 6–10 2022
```

