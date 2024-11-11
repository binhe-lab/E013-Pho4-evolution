Use PADDLE (Sanborn et al 2021 eLife) to predict activation potential in ScPho4 and CgPho4.

The PADDLE scripts and models were downloaded from <https://github.com/asanborn/PADDLE> on 2024-03-05

To install `tensorflow` on my M1 mac, I used the following lines

```unix
conda create -n E013 python=3.10.0 tensorflow
conda deactivate; conda activate E013
mamba install -c conda-forge jupyterlab
pip install matplotlib pandas
```
