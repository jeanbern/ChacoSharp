## ChacoSharp
A C# port of the [Chaco](https://github.com/gsjaardema/seacas/tree/master/packages/seacas/libraries/chaco) graph partitioning library.  
This project may be incomplete or even broken. I've only used it for spectral sequencing.

### Documentation
* There is a [pdf guide](https://github.com/jeanbern/ChacoSharp/blob/master/ChacoUserGuide2.0.pdf) available detailing some of the options and usage flags available.  
* The aforementioned configurations are defined in [StaticConstants.cs](https://github.com/jeanbern/ChacoSharp/blob/master/StaticConstants.cs), often using descriptive enum types or bool instead of cryptic integers. Many of these are *const* and will require re-compilation of both the project and its consumers for changes to apply.  

### Credits
* The original code is copyright (c) 2005-2017 National Technology & Engineering Solutions of Sandia, LLC (NTESS). It can be found [here](https://github.com/gsjaardema/seacas/tree/master/packages/seacas/libraries/chaco)
* *ChacoSequenceJPetit.cs* is a replacement for [sequence.c](https://github.com/gsjaardema/seacas/blob/master/packages/seacas/libraries/chaco/misc/sequence.c) as found in the [work of Jordi Petit](https://www.cs.upc.edu/~jpetit/MinLA/Experiments/).
* Porting to C# done by [Jean-Bernard Pellerin](https://jbp.dev) - [github](https://github.com/jeanbern)
* In the absence of a copyright claim or license, I've used Sandia's [Guide to Preparing SAND Reports and Other Communication Products](https://github.com/jeanbern/ChacoSharp/blob/master/SandiaReportGuide.pdf) to interprete *"Unlimited Release"* in [The Chaco User's Guide Version 2.0](https://github.com/jeanbern/ChacoSharp/blob/master/ChacoUserGuide2.0.pdf)