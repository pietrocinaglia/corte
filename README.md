# tGeCoNet: a framework for constructing (t)emporal (Ge)ne (Co)-expression (Net)works

**tGeCoNet** is a python software package for constructing (real-world) temporal networks from age-related gene expression data (i.e., temporal gene co-expression networks).


## Additional Information

Age-related gene expression data is retrieved by Genotype-Tissue Expression (GTEx) API V2 (https://gtexportal.org)

Time points were based on the following age ranges: ['20-29', '30-39', '40-49', '50-59', '60-69', '70-79']


## Table of Contents

- [Installation](#installation)
- [Usage](#usage)
    - [Constructing a Temporal Network from age-related gene expression data](#constructing)
    - [Plotting](#plotting)
- [Resuts](#results)
- [License](#license)


## Installation <a name="installation"></a>

To begin using tGeCoNet you have to follow these steps:

1. Install dependecies

```
pip install -r requirements.txt
```

2. Import tGeCoNet in your code:

- 2.1. Copy "tgeconet" folder into your workspace.

- 2.2. Import tGeCoNet in your code:

```
from tgeconet import TGECONET
```

3. Create an instance of **tGeCoNet**:

```
tgeconet = TGECONET(genes_of_interest:list, tissues_of_interest:list, threshold:float=0.05, verbose:bool=False):
```

4. Use the following methods:

- construct_temporal_network() -> list
- plot(with_labels:bool=True, savefig:str=None)

Additionally, the following method for testing was implemented:

- test(plot:bool=True, savefig:str=None)

## Usage <a name="usage"></a>

**tGeCoNet** allows for constructing (real-world) temporal networks from age-related gene expression data

You can choose the method that best suits your data and use case.


# Constructing a Temporal Network from age-related gene expression data <a name="constructing"></a>

```
tgeconet = TGECONET(genes_of_interest:list, tissues_of_interest:list, threshold:float=0.05, verbose:bool=False):
temporal_network = tgeconet.construct_temporal_network(self)
```

# Plotting the temporal network <a name="plotting"></a>

'output_path' is the path where files will be saved.
Note that if 'savefig' is None, temporal network will only be shown.

```
tgeconet.plot(with_labels:True, savefig=output_path)
```


## Results <a name="results"></a>

The 'results' folder contains the results from tests conducted in our own experimentation and discussed in the academic article related to **tGeCoNet**


### License <a name="license"></a>

MIT License - Copyright (c) Pietro Cinaglia

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.