# spatial-epi
A collection of bits and pieces for thinking through how to do spatial epidemic spread modelling.

In particular a series of Netlogo models, as follows. Two (very) toy models exploring self isolation 'bubbles'
+ [`bubbles.nlogo`](http://southosullivan.com/misc/bubbles.html)
+ [`nested-bubbles.nlogo`](http://southosullivan.com/misc/nested-bubbles.html)

and more substantively, a series of connected localised SEIR models. These have been 'frozen' for reference purposes and because 'releases' aren't really appropriate to this project. Brief details as follows, with links to web version where available:
+ [`distributed-seir.nlogo`](http://southosullivan.com/misc/distributed-seir.html) has random connections among a set of equal-sized locales
+ [`distributed-seir-02.nlogo`](http://southosullivan.com/misc/distributed-seir-02.html) has more spatially coherent connections among the same
+ `distributed-seir-03.nlogo` is an attempt to enable the model to read in real health management zones. It works, but requires some idiosyncratic code for the file reading.

All three of the above have excessive mortality, which has been corrected in later models. The later models are more worth spending time with.
+ `distributed-seir-04.nlogo` reverts back to the main sequence and allows for locales to vary in size and variance, under paramterisable control
+ `distributed-seir-05.nlogo` fixed lockdown levels and testing added
+ [`distributed-seir-06.nlogo`](http://southosullivan.com/misc/distributed-seir-06.html) as previous but with correction to locale sizes to match total population more closely
+ [`distributed-seir-07.nlogo`](http://southosullivan.com/misc/distributed-seir-07.html) as previous but with automatic control of lockdown levels according to a variety of strategies
+ [`distributed-seir-08.nlogo`](http://southosullivan.com/misc/distributed-seir-08-web.html) adds logging of all locales data at every time step (the linked web version has no file logging)

And now a version that can also be initialised with NZ DHB data:
+ [`nz-dhb-seir-08.nlogo`](http://southosullivan.com/misc/nz-dhb-seir-08-web.html) adds logging reading of spatial data from input GUI elements (done this way to permit same in web version)

COVID-19 spread parameters based on values used in [this work](https://cpb-ap-se2.wpmucdn.com/blogs.auckland.ac.nz/dist/d/75/files/2017/01/Supression-and-Mitigation-Strategies-New-Zealand-TPM-1.pdf), although results unlikely to match exactly given entirely different platform used (and the rapidly evolving situation).

You can make a web version of any of these by uploading to [http://netlogoweb.org/](http://netlogoweb.org/launch#Load)
