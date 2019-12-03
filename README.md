## ChacoSharp
A C# port of the [Chaco](https://github.com/gsjaardema/seacas/tree/master/packages/seacas/libraries/chaco) graph partitioning library.  
This project may be incomplete or even broken. I've only used it for spectral sequencing.

### Documentation
* There is a [pdf guide](https://gsjaardema.github.io/seacas/chaco.pdf) available detailing some of the configuration options.  
* The flags mentioned can be set in *StaticConstants.cs*, often using descriptive enum types or bool instead of cryptic integers.
* See [placeholder link here](github.com/jeanbern) for example usage.

### Credits
* The original code is copyright (c) 2005-2017 National Technology & Engineering Solutions of Sandia, LLC (NTESS). It can be found [here](https://github.com/gsjaardema/seacas/tree/master/packages/seacas/libraries/chaco)
* *ChacoSequenceJPetit.cs* is a replacement for *[sequence.c](https://github.com/gsjaardema/seacas/blob/master/packages/seacas/libraries/chaco/misc/sequence.c)* as found in the [work of Jordi Petit](https://www.cs.upc.edu/~jpetit/MinLA/Experiments/).
* Porting to C# done by [Jean-Bernard Pellerin](jbp.dev) - [github](github.com/jeanbern)