# VEmodel
## Repository for vibronic models construction.
<div align="center">

![VEmodel_scheme](https://github.com/user-attachments/assets/65604285-ce71-4c55-b9ea-ffd212dab3e2)

<h3>

</h3>

</div>

---

### Running the program (Toby_automate folder)

```py
python3 dist_allmodes_pm.py nh3_ccd_mp2_c3v_gh.out
```

__nlogn__: change GAMESS submission script path to _/home/bjb2chen/LOCAL/subgam.diab_

__ComputeCanada(Graham)__: change path to _/home/bjb2chen/LOCAL/subgam.mrsftd_ (or contact tzeng@yorku.ca)

Video recording on how to run the diabatization code (assuming you have GAMESS with gmcpt module already): 

- [x] C2v H2O diabatization - [https://www.youtube.com/watch?v=ILLrUSzR3YQ&ab_channel=BennyChen]

- [x] D3h NH3 diabatization - [https://www.youtube.com/watch?v=zeKTlq3VHho&ab_channel=BennyChen]

- [x] D3h CoF3 diabatization - [https://www.youtube.com/watch?v=7tNUAFYKCQE&ab_channel=BennyChen]

- [x] [Project Poster Reference](https://github.com/bjb2chen/vmodels/files/10171706/SCP2022_bjc_20685630_White.pdf)

![nice_protocol](https://github.com/bjb2chen/VEmodel/assets/51763900/8ea2d42e-e0ec-4aac-b2f2-47cb1d31348e)

## Vibronic Spectra - Test Case Molecules

### $H_2O^+$

![H2O_waterfall](https://github.com/bjb2chen/VEmodel/assets/51763900/a793d492-7580-487f-81a2-4ff793b3f7c2)

### $NH_3^+$

![NH3_waterfall](https://github.com/bjb2chen/VEmodel/assets/51763900/9d201627-a57e-483c-9850-f920759051c5)

### $PH_3^+$

![Jun_12_PH3_composite](https://github.com/bjb2chen/VEmodel/assets/51763900/12dd6e52-24f2-4f1d-861d-842f0caad53f)

### $CH_2O^+$

![Jun_13_CH2O_composite](https://github.com/bjb2chen/VEmodel/assets/51763900/af4d06a3-d6da-4917-97c7-8f07cbab10a5)

## Vibronic Spectra - Transition Metal Complexes

### $D_{3h}$ $RhF_3$

(tf = 2000fs, tau = 2000, iexp = 1)
![RhF3_spectrum](https://github.com/bjb2chen/VEmodel/assets/51763900/97234174-5bd0-4f83-a691-374efb1cabd3)

![Jul_08_RhF3_constantsH0S0](https://github.com/bjb2chen/VEmodel/assets/51763900/fe67e2b4-c3b1-4ce4-8496-11782e9b101c)
![Jul_08_RhF3_Linears](https://github.com/bjb2chen/VEmodel/assets/51763900/e4a96feb-a74e-40cf-90bf-65432419748a)
![Jul_08_RhF3_Quadratics](https://github.com/bjb2chen/VEmodel/assets/51763900/3f78443d-941e-4064-9a6d-84e87573a579)

### $D_{3h}$ $CoF_3$

![Jul_08_CoF3_constants](https://github.com/bjb2chen/VEmodel/assets/51763900/fa47ffd0-13e0-443e-8d23-93aad02f73b0)
![Jul_08_CoF3_Linears](https://github.com/bjb2chen/VEmodel/assets/51763900/248b6df7-2046-43ac-8f8f-cfaa86640c8b)
![Jul_08_CoF3_Quadratics](https://github.com/bjb2chen/VEmodel/assets/51763900/0107f704-7b9d-486c-be1a-c0e51fe5b54f)

### $Fe(CO)_5$

![Jun_25_FeCO_waterfall_composite_collated](https://github.com/bjb2chen/VEmodel/assets/51763900/83e94785-24c1-4635-88dc-fe145d9fc277)

<img width="1013" alt="waterfall_FeCO" src="https://github.com/bjb2chen/VEmodel/assets/51763900/09098d89-2252-4fe2-80d8-00c8af03f0d2">

![FeCO_PBF5_200fs_coupledmodes](https://github.com/bjb2chen/VEmodel/assets/51763900/a3f18842-0c68-46ff-b59f-cd49ec27553a)

![CoF4](https://github.com/bjb2chen/vmodels/assets/51763900/eb5d7752-d0d4-4151-9af5-d399e079bf3a)



