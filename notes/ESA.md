# 一维电扫阵列

## 方向图公式

- ![image-20200613104318957](C:\Users\sfpotato\AppData\Roaming\Typora\typora-user-images\image-20200613104318957.png)
- $EP=cos^{\frac{EF}{2}} \theta$，其中$EF$表示阵元因子，$EP$表示阵元方向图

![image-20200613104551942](C:\Users\sfpotato\AppData\Roaming\Typora\typora-user-images\image-20200613104551942.png)

## 基本参数

- 波束宽度:$\theta_{bw} = \frac{k\lambda}{Mdcos\theta_0}$，其中*k*表示波束宽度因子，当$k=0.886$时且为均匀口径照射，波束宽度为3dB带宽

- 瞬时带宽：瞬时带宽表示损耗能够在接受范围内的频率变化范围

  $IBW=\frac{kc}{Lsin\theta_0}$

- 栅瓣表示主瓣以外的气压方向有规律地形成与主播书类似地辐射波束，栅瓣的位置与频率和阵元间距有关

  ![image-20200613105606723](C:\Users\sfpotato\AppData\Roaming\Typora\typora-user-images\image-20200613105606723.png)

- 对于均匀的幅值分布，阵元方向图不发生扫描，阵因子发生扫描，导致增益下降，峰值位置产生偏差，但是均匀加权的方式是能量最高的![image-20200613105809290](C:\Users\sfpotato\AppData\Roaming\Typora\typora-user-images\image-20200613105809290.png)

# 二维电扫阵列

## 方向图公式

- ![image-20200613115843143](C:\Users\sfpotato\AppData\Roaming\Typora\typora-user-images\image-20200613115843143.png)

  坐标表示：

  $x_m = (m - \frac{M+1}{2})d_x, m = 1, ..., M; x_y = (n - \frac{N+1}{2})d_y, n = 1, ..., N$

  ![image-20200613120242003](C:\Users\sfpotato\AppData\Roaming\Typora\typora-user-images\image-20200613120242003.png)

  ![image-20200613120312504](C:\Users\sfpotato\AppData\Roaming\Typora\typora-user-images\image-20200613120312504.png)

  

## ESA的空间坐标定义

- 天线坐标系

  ![image-20200613120449873](C:\Users\sfpotato\AppData\Roaming\Typora\typora-user-images\image-20200613120449873.png)

- 雷达坐标系

  ![image-20200613120532782](C:\Users\sfpotato\AppData\Roaming\Typora\typora-user-images\image-20200613120532782.png)

- 天线圆锥坐标系

  ![image-20200613120608195](C:\Users\sfpotato\AppData\Roaming\Typora\typora-user-images\image-20200613120608195.png)

- 三个坐标系之间的转换

  ![image-20200613120644174](C:\Users\sfpotato\AppData\Roaming\Typora\typora-user-images\image-20200613120644174.png)

  其中$arctan2(y, x) = arctan(y/x)$

- 正弦空间表示法

  ![image-20200613122127456](C:\Users\sfpotato\AppData\Roaming\Typora\typora-user-images\image-20200613122127456.png)

  ![image-20200613122211231](C:\Users\sfpotato\AppData\Roaming\Typora\typora-user-images\image-20200613122211231.png)

- ![image-20200613122314824](C:\Users\sfpotato\AppData\Roaming\Typora\typora-user-images\image-20200613122314824.png)

  在这里波束宽度与扫描角无关，扫描波束的峰值为$sin\theta_z$

## 阵元网格

- 栅瓣扫描盲区，只有当主波束处于栅瓣扫描盲区时，栅瓣才不会出现在真实空间中，即以主波束为中心的圆内

  ![image-20200613134853177](C:\Users\sfpotato\AppData\Roaming\Typora\typora-user-images\image-20200613134853177.png)

- ![image-20200613135947750](C:\Users\sfpotato\AppData\Roaming\Typora\typora-user-images\image-20200613135947750.png)

- 阵元的三角形分布

  ![image-20200613140250273](C:\Users\sfpotato\AppData\Roaming\Typora\typora-user-images\image-20200613140250273.png)

- ![image-20200613140328863](C:\Users\sfpotato\AppData\Roaming\Typora\typora-user-images\image-20200613140328863.png)

## 二维阵列方向图分析

- 考虑ESA的误差影响：①相关误差：移相器和衰减器的量化误差②故障元件以及制造偏差等

- 量化误差引起的副瓣电平峰值和平均值计算公式：

  $Average\ SLL = \frac{\pi ^2}{3n_{elem} 2^{2N}}; Max \ SLL = \frac{1}{2^{2N}}$，其中$N$表示量化位数。

- 俯仰角、横滚角和偏航角引起的坐标变化可以用旋转矩阵来表示

  ![image-20200614093225010](C:\Users\sfpotato\AppData\Roaming\Typora\typora-user-images\image-20200614093225010.png)

  <img src="C:\Users\sfpotato\AppData\Roaming\Typora\typora-user-images\image-20200614093401725.png" alt="image-20200614093401725" style="zoom:80%;" />

  

- 

# 子阵列波束形成

## 概念

- ![image-20200624103353525](C:\Users\sfpotato\AppData\Roaming\Typora\typora-user-images\image-20200624103353525.png)
- ![image-20200624103409158](C:\Users\sfpotato\AppData\Roaming\Typora\typora-user-images\image-20200624103409158.png)
- 只需要在各级子阵后接入移相器就行了