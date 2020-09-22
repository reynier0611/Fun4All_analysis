On the Berkeley EIC meeting on Sept 15, 2020, I presented a graph of DCA_z vs. momentum, in which I had three different detector settings:
1) original geometry by Ernst Sichtermann and Yue Shi Lai (2 vertexing layers, 0.3% X/X0),
2) updated geometry to account for new projected EIC beampipe (just like in the previous configuration, but with larger radii), and
3) new geometry (3 vertexing layers, 0.05% X/X0)

It was correctly pointed out that comparing 1 or 2 vs. 3 was not straightforward, since in configuration 3 two parameters
(# layers and stave material budget) had changed.
In this directory, I study a simple geometry in which I created the barrel layers from the All-Si tracker geometry in a simplified way.
Each layer is a perfectly smooth cylinder (rather than a cylinder made of staves). There are three vertexing layers of 0.05% X/X0, and the
remaining tracking layers are of 0.5% X/X0. The beryllium portion of the beampipe is also created.
The code can be found [here](https://github.com/reynier0611/g4lblvtx/tree/master/macros/auxiliary_studies/vertexing).
In this code, I can turn on or off different vertexing layers. Which layers are on or off are represented by a boolean. For instance:
* 111 means that all three vertexing layers are turned on
* 110 means that the inner-most and middle vertexing layers are turned on, while the outer-most vertexing layer is off.
The codes in this directory are meant to plot the output of the files produced with this code.
