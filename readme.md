## 実行の仕方
1. `cd fuzzy-broccoli`
2. `make`
3. `bin/main`

各種ファイルは`fuzzy-broccoli`からの相対パスで指定されているため,`cd bin`,`./main`としてしまうとエラーが発生する.

## コンパイル途中に生成されるファイルの削除
1. `make clean`

## 考慮している現象
水,空気の二相流と塩輸送の連成解析を考慮している.

## 考慮している問題
現時点では、以下の問題を考慮している。詳しくはReferenceを参照。
圧力120bar、温度 45◦C、塩の質量濃度15%の環境で、厚さ100mの均質で等方性のある帯水層を貫通している井戸がある。解析領域は、解析期間中に影響を与えない程度に十分広く(100km)設定する。二酸化炭素は100kg/sの一定速度で均一に圧入される。

## 計算方法
有限差分法,１次風上法,陰解法を用いている.
非線形計算にはNewton-Rapson法を用いている.
偏微分係数の計算には複素数階微分法を用いている。

## Reference
Karsten Pruess. ECO2N: A TOUGH2 fluid property module for mixtures of water, NaCl, and CO2. Lawrence Berkeley National Laboratory Berkeley, CA, 2005.を参照
