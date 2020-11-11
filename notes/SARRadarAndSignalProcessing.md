# 合成孔径雷达介绍

#### 参考书籍

- 《合成孔径雷达成像-算法与实现》

**课程安排：四次课程设计**

-------

### Notes

- 天波超视距雷达使用3MHz-30MHz的短波波段
- 点迹处理：将同一目标的点匹配
- 航迹处理：将不同*PRT*的目标形成的点处理形成目标的运动路径
- “大气窗口”，选择大气透射率较高的波段作为雷达的波段

![image-20200506144722759](C:\Users\sfpotato\AppData\Roaming\Typora\typora-user-images\image-20200506144722759.png)

- 垂直极化不是指与地面垂直的极化方向，而是与水平极化和*k*构成的平面右手法则的方向
- 孔径越大，波束越窄，分辨率越高
- 积分旁瓣比，峰值旁瓣比（主瓣与第一旁瓣（旁瓣中峰值最大旁瓣而不是离主瓣最近的旁瓣）的分贝值之差）
- 采用脉冲压缩可以同时提高距离分辨率和占空比：$dR = \frac{c}{2B}$
- 定位精度与分辨率，如果场景中有多个目标定位精度等同于分辨率，否则定位精度可通过多个测量数据超过分辨率

### SAR的成像基本原理

- ![image-20200506153358528](C:\Users\sfpotato\AppData\Roaming\Typora\typora-user-images\image-20200506153358528.png)

- 成像特点：叠影（距离接收机相同距离的点会叠在一起），前景压缩，阴影。

- 合成孔径的原理图（D为真实孔径，SAR的分辨率为$\frac{D}{2}$）

  ![image-20200506153907649](C:\Users\sfpotato\AppData\Roaming\Typora\typora-user-images\image-20200506153907649.png)

- SAR的工作模式

  ![image-20200506154348196](C:\Users\sfpotato\AppData\Roaming\Typora\typora-user-images\image-20200506154348196.png)

- 运动补偿方法：①自聚焦方法②基于IMU/GPS的运动补偿

-----

# 信号基础

### 信号的描述与分类

- 从一个脉冲看是$f(t)$一维信号，从多个脉冲看是多维信号

- 利用矩形脉冲的思想实现傅里叶变换（从理论到实际）

  ![image-20200509143618177](C:\Users\sfpotato\AppData\Roaming\Typora\typora-user-images\image-20200509143618177.png)

  坐标转换

  ![image-20200509144033076](C:\Users\sfpotato\AppData\Roaming\Typora\typora-user-images\image-20200509144033076.png)

### 信号与系统的概念

- 时域分析法优点是精确，频域分析法优点是快

  #### 卷积性质

  - ![image-20200509144811188](C:\Users\sfpotato\AppData\Roaming\Typora\typora-user-images\image-20200509144811188.png)
  - 长为M, N的序列卷积后的结果长度为$M + N - 1$

  #### 傅里叶变换性质

  - $FT[F(w-w_0)] = f(t) e^{jw_0t}$
  - $FT[f(t-t_0)] = F(w) e^{-jwt_0}$
  - ![image-20200509150512936](C:\Users\sfpotato\AppData\Roaming\Typora\typora-user-images\image-20200509150512936.png)
  - $f(nT_s) = \frac{1}{T_s} \sum_{k=0}^{N-1} F(kw_1)e^{jknw_1T_s}$

# 线性调频信号分析与脉冲压缩

### 信号抽样

- ![image-20200509152357654](C:\Users\sfpotato\AppData\Roaming\Typora\typora-user-images\image-20200509152357654.png)

  数字信号的频谱由于采样保持$u(t)$ 即$FT[f(t)\otimes u(t)] = F(w)*sinc(w)$导致频谱畸变

### 线性调频信号

- ![image-20200513120939341](C:\Users\sfpotato\AppData\Roaming\Typora\typora-user-images\image-20200513120939341.png)

- 时域表达式$s(t) = rect(\frac{t}{T}) exp \left[j \pi K t^2 \right]$

  频域表达式$S(f) = \int _{-\infty}^{\infty}rect(\frac{t}{T}) exp \left[j \pi K t^2 \right]exp\left[-j2\pi ft\right]dt$

  难以直接利用傅里叶变换直接推导，可利用驻定相位原理得到近似表达式

- 驻定相位原理 ：对于$g(t) = w(t) exp[j \phi (t)]$的积分，信号在相位驻留点附近是缓变的，而在其他时间点上事捷变的，驻留点附近对积分起主要贡献(只要$TBP$够大，驻定相位原理是相当准确的)

  ![image-20200513120202004](C:\Users\sfpotato\AppData\Roaming\Typora\typora-user-images\image-20200513120202004.png)

  ![image-20200513120615855](C:\Users\sfpotato\AppData\Roaming\Typora\typora-user-images\image-20200513120615855.png)

- 线性调频信号与采样间隔的关系，过采样率保持再1.2以上才可以保持不失真

- 匹配滤波器的物理特点：幅频特性与信号一致；相频特性与信号相反，使得信号同相叠加

- 时域翻转求共轭对应频域直接求共轭

- 匹配滤波器：$h(t) = ks^*(t_m - t)$，冲激响应是发射信号的时域翻转和共轭

  $H(jw) = k[S(jw)e^{jwt_m}]^*$

### 脉冲压缩

- 基本思想：补偿相位

  ![image-20200509155754092](C:\Users\sfpotato\AppData\Roaming\Typora\typora-user-images\image-20200509155754092.png)

- 时域实现：
  - ![image-20200516101220386](C:\Users\sfpotato\AppData\Roaming\Typora\typora-user-images\image-20200516101220386.png)

  当$TBP>100$时，$S_{out} (t) \approx Tsinc[KT(t - t_0)]$

  LMF经过匹配滤波后

  ![image-20200516104427092](C:\Users\sfpotato\AppData\Roaming\Typora\typora-user-images\image-20200516104427092.png)

  ![image-20200516104730241](C:\Users\sfpotato\AppData\Roaming\Typora\typora-user-images\image-20200516104730241.png)

  时间量纲下的3dB分辨率： $\rho = \frac{0.886}{|K|T} \approx \frac{1}{|K|T}$

  压缩比为时间带宽积：$r = \frac{\rho}{\rho ^{\prime}} = \frac{T}{1/(|K|T)} = |K|T^2$

  SNR提高倍数为压缩比

  脉冲压缩后的时间分辨率为信号带宽的倒数

  LMF压缩后的信号为实数。

  ![image-20200516111400087](C:\Users\sfpotato\AppData\Roaming\Typora\typora-user-images\image-20200516111400087.png)

  相位为$\pi,0$

  -----

- 频域实现

  线性调频信号在带宽内均匀扫频，具有平坦的频谱，频域中的脉冲压缩本质就是将信号频谱与含有二次共轭相位 的频域滤波器相乘，以得到具有线性相位的平坦频谱

  基带信号

  ![image-20200516112303620](C:\Users\sfpotato\AppData\Roaming\Typora\typora-user-images\image-20200516112303620.png)

  非基带信号

  ![image-20200516113139820](C:\Users\sfpotato\AppData\Roaming\Typora\typora-user-images\image-20200516113139820.png)

  频域匹配滤波器生成方式

  ![image-20200516113451459](C:\Users\sfpotato\AppData\Roaming\Typora\typora-user-images\image-20200516113451459.png)

- 峰值旁瓣比与积分旁瓣比

  ![image-20200516114607718](C:\Users\sfpotato\AppData\Roaming\Typora\typora-user-images\image-20200516114607718.png)

  如何降低PSLR：加权（加窗），使得PSLR降低，但同时引起波束宽度增加，SNR降低

- **TBP**：时间带宽积

- 如何插值：在频域两边补零，然后傅里叶变换

### Q&A

- Q信号带宽

  A例如由三个子基带构成的信号，如果需要3 dB带宽为100 M，则需要每个子基带的1 dB带宽为100 M

# SAR成像

### SAR几何关系

- 波束宽度分为方位向和距离向

  斜视角表示SAR天线的朝向；零多普勒面表示多普勒速度为0的位置，位于这个平面内的目标与飞机的径向速度为0

  ![image-20200511145319798](C:\Users\sfpotato\AppData\Roaming\Typora\typora-user-images\image-20200511145319798.png)

- 距离等式的双曲线模型：$R^2 (\eta) = R_0 ^2 + V_r ^ 2 \eta ^2$

  其中$R^2 (\eta)$表示斜距，$R_0$表示零多普勒面的斜距，$V_r$表示飞机飞行速度

  利用泰勒公式进行近似

  ![image-20200516115536916](C:\Users\sfpotato\AppData\Roaming\Typora\typora-user-images\image-20200516115536916.png)

  其中$\theta _c$表示零多普勒面与斜视的夹角

  小斜视角可以近似为 $R(\eta) \approx R_0 + \frac{V^{2}\eta ^2}{2R_0}$

### 距离向信号分析

- 接收信号的带宽为 $B=|K_r|T_r$

  距离向的分辨率（与脉冲宽度有关）$\rho _r = \frac{c}{2|K_r| T_r}$
  
- 距离向的采样率要保证能够准确刻画波形

###  方位向信号分析

- 方位向的波束扫描

  ![image-20200516122146334](C:\Users\sfpotato\AppData\Roaming\Typora\typora-user-images\image-20200516122146334.png)

  ![image-20200516122717892](C:\Users\sfpotato\AppData\Roaming\Typora\typora-user-images\image-20200516122717892.png)

  由四部分构成：发射信号的加权函数$w_r$，天线方向图$w_a$（双程，类似于$sinc^2$），载频的多普勒频移，LFM的多普勒频移

- 方位向的多普勒参数

  方位向的多普勒速度

  ![image-20200516144444589](C:\Users\sfpotato\AppData\Roaming\Typora\typora-user-images\image-20200516144444589.png)

  ![image-20200516145436357](C:\Users\sfpotato\AppData\Roaming\Typora\typora-user-images\image-20200516145436357.png)

  $\theta _{bw}$表示波束宽度，多普勒带框可以理解为距离向的两端的多普勒频移的差值

  ![image-20200516145814427](C:\Users\sfpotato\AppData\Roaming\Typora\typora-user-images\image-20200516145814427.png)

  调频率$K_a$是负数，表明多普勒频率由负到正

- 方位向分辨率

  ![image-20200516150655537](C:\Users\sfpotato\AppData\Roaming\Typora\typora-user-images\image-20200516150655537.png)
  
- 合成孔径:$D_s = R_0\lambda / D$

### PRF的选择

- 测绘带与方位向分辨率对PRF的要求存在悖论：方位向分辨率越高要求多普勒带宽越大，那么PRF就得高即$f_c$高，那么PRT就得小，那么测绘带就得小
- PRF其实就是方位向的采样率
- ![image-20200517175033380](C:\Users\sfpotato\AppData\Roaming\Typora\typora-user-images\image-20200517175033380.png)
- $PRT = T_r + T_r + 2(R_2 - R_1)/c$表示为反射脉冲的宽度 + 最后一个接收回波信号的脉冲宽度 + 最近斜距与最远斜距之差

-----

# SAR成像算法

- 对于“成像”可以有两种理解：
  - 将二维 （方位向、距离向）的弥散的回波信号进行距离向和方位向的二维压缩得到点信号，即回波脉冲压缩
  - 对回波进行反卷积从而获得目标的点信息

### RD算法

- Range-Doppler算法：距离徙动矫正 > 距离向压缩 > 方位向压缩

  ![image-20200518113600482](C:\Users\sfpotato\AppData\Roaming\Typora\typora-user-images\image-20200518113600482.png)

  由于方位向的多普勒频率是线性调频的，所以可以将不同点目标的回波信号进行距离向的压缩

  距离徙动矫正的两种表达方式

  ![image-20200518120154669](C:\Users\sfpotato\AppData\Roaming\Typora\typora-user-images\image-20200518120154669.png)

  ![image-20200518120237084](C:\Users\sfpotato\AppData\Roaming\Typora\typora-user-images\image-20200518120237084.png)

  其中$f_{dc}$表示多普勒中心频率，$f_{dr}$表示多普勒调频率。线性项表示距离走动，二次项表示距离

- 距离向脉冲压缩后在距离多普勒域内常常出现卷绕的现象

  ![image-20200518141532205](C:\Users\sfpotato\AppData\Roaming\Typora\typora-user-images\image-20200518141532205.png)

  ==原因在于多普勒频率中心$f_{dc}$不为0==

- ![image-20200518151702491](C:\Users\sfpotato\AppData\Roaming\Typora\typora-user-images\image-20200518151702491.png)

- 估计多普勒中心解决多普勒模糊问题

# 基于运动补数据运动补偿

- ![image-20200629114841773](C:\Users\sfpotato\AppData\Roaming\Typora\typora-user-images\image-20200629114841773.png)

  ### 匀速误差

- ![image-20200629114907887](C:\Users\sfpotato\AppData\Roaming\Typora\typora-user-images\image-20200629114907887.png)

  其中$r_1(\eta) - r(\eta)$省去了$\Delta V_A$的二次项

  ### 匀加速误差

  ![image-20200629115308795](C:\Users\sfpotato\AppData\Roaming\Typora\typora-user-images\image-20200629115308795.png)

  ### 正弦速度误差

  ![image-20200629115438044](C:\Users\sfpotato\AppData\Roaming\Typora\typora-user-images\image-20200629115438044.png)

  当存在与飞行方向垂直的横向速度分量时，产生雷达视线方向运动损失，导致基带线性调频信号变成非基带线性调频![image-20200629115716852](C:\Users\sfpotato\AppData\Roaming\Typora\typora-user-images\image-20200629115716852.png)

  ![image-20200629120623034](C:\Users\sfpotato\AppData\Roaming\Typora\typora-user-images\image-20200629120623034.png)

  - 距离重采样：频域乘上线性相位项
  - 方位重采样：sinc插值

- ![image-20200629121305852](C:\Users\sfpotato\AppData\Roaming\Typora\typora-user-images\image-20200629121305852.png)

- ![image-20200629113730305](C:\Users\sfpotato\AppData\Roaming\Typora\typora-user-images\image-20200629113730305.png)

  第一步针对参考斜距进行一阶运动补偿

  ![image-20200629113809795](C:\Users\sfpotato\AppData\Roaming\Typora\typora-user-images\image-20200629113809795.png)

  $z_i, x_i$是飞机的坐标

  #### 一阶运动补偿

  ![image-20200629113842455](C:\Users\sfpotato\AppData\Roaming\Typora\typora-user-images\image-20200629113842455.png)

  #### 距离重采样

  在距离向频域乘上距离向的相位项（只对中心测绘带的误差进行补偿）

  ![image-20200629114002071](C:\Users\sfpotato\AppData\Roaming\Typora\typora-user-images\image-20200629114002071.png)

  #### 方位向重采样

  一点一点算sinc插值

  ![image-20200629140107962](C:\Users\sfpotato\AppData\Roaming\Typora\typora-user-images\image-20200629140107962.png)

  ![image-20200629223043980](C:\Users\sfpotato\AppData\Roaming\Typora\typora-user-images\image-20200629223043980.png)

  

  

# 基于回波数据的参数估计与运动补偿

### 多普勒中心频率估计

- 一阶相位误差

- 当波束方向图主瓣没有对准零多普勒中心时，导致模糊比增高![](C:\Users\sfpotato\AppData\Roaming\Typora\typora-user-images\image-20200601144152201.png)
- 方位模糊比：模糊信号能量之和与主瓣能量之和的比值
- 将每个多普勒频谱点能量用指数概率密度函数进行建模，采用最大似然估计，得到最优滤波器

![image-20200601150506559](C:\Users\sfpotato\AppData\Roaming\Typora\typora-user-images\image-20200601150506559.png)

- 方向图越陡多普勒中心估计越精确
- ![image-20200601151802679](C:\Users\sfpotato\AppData\Roaming\Typora\typora-user-images\image-20200601151802679.png)
- 圆周卷积一般用傅里叶域相乘实现，注意傅里叶变换不要补零

### 多普勒调频率估计

- 二次相位误差

- ![image-20200601155455456](C:\Users\sfpotato\AppData\Roaming\Typora\typora-user-images\image-20200601155455456.png)

- ![image-20200601155440610](C:\Users\sfpotato\AppData\Roaming\Typora\typora-user-images\image-20200601155440610.png)

### 相位梯度自聚焦算法

- 高次相位误差

- ![image-20200601160159766](C:\Users\sfpotato\AppData\Roaming\Typora\typora-user-images\image-20200601160159766.png)
- ![image-20200601160520819](C:\Users\sfpotato\AppData\Roaming\Typora\typora-user-images\image-20200601160520819.png)

