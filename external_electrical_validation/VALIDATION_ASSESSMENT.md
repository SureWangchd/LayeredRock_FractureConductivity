# 外部观测验证判定

**结论：B，部分验证。已获得真实、可复算的实验—预测比较；未达到核心主张的完整独立实验验证。**

|判定项目|实际完成情况|
|---|---|
|真实观测|Han2022表1的60条电学/孔隙率数据已核实；固定5 MPa使用六个不同岩样|
|校准与目标分离|只用I-0/I-90标定；I-45与三个F样品电导率全部留出|
|实际科学计算|贴体三维有限元；有限椭球单胞加有限圆柱端面电极算子；四档网格|
|总电导率|I-45误差−1.72%；三个F样品原论文闭合误差1.38%–2.03%，基质域闭合2.51%–3.11%|
|方向排序|两种闭合均得到F-0 < F-45 < F-90；样品孔隙率不同，排序不能单独归因于裂缝取向|
|裂缝增量|从跨样品背景转移得到，未直接测量；点值误差大于总量误差，且依赖孔隙率闭合|
|几何独立性|制造尺寸、密度、取向来自实验；压力下开度与逐裂缝坐标缺失，使用固定规则阵列与情景|
|材料独立性|F目标未参与标定；盐水电导率跨研究转移，已测试独立量级情景；同批实测值未提供|
|测量物理|原观测为1 kHz；端面条件已对应；普通介电项很小，但相位/电极极化未获独立检验|
|原稿核心命题|背景匹配后的结构重排、外围干预及模态预测没有相应直接实验|

## 能支持的文字

“Using intact-specimen electrical measurements for background calibration, the finite-element model reproduces the measured ordering and approximates the total conductivities of three fractured artificial sandstone specimens at 5 MPa. The comparison uses the reported 1 kHz conductive response. Its accuracy depends on the transfer of background porosity and fracture geometry.”

模型数值误差、实验误差与模型情景分别报告。两种孔隙率闭合的预测均保留，未按目标吻合程度选择一个作为真实重建。

## 本数据不能支持的文字

- “The matched-background A/B effect has been experimentally validated.”
- “All fracture-induced increments are predicted within 2%.”
- “The 60 observations provide 60 independent validation samples.”
- “The measured 5 MPa crack apertures were reconstructed.”
- “The 1 kHz data prove a frequency-independent DC plateau.”
- “The exterior modal correction has been independently validated by laboratory observations.”

## 三项清楚的证据层级

1. **已经完成：** 公开人工砂岩总电导率的外部定量比较；真实求解，未拟合F目标。
2. **当前不能唯一识别：** 三方向裂缝净增量。需要同岩样前后对照，或独立约束的基质背景、压力下裂缝几何及同批盐水/测量频谱。
3. **不可替代的主命题实验：** 同背景等效条件、同一有限裂缝或固定含裂缝核心，改变内部排列/外围结构，测得不同新增电导。再用未参与标定的状态检验外围响应模型。

因此，对“是否真正补上JGRSE的观测验证缺口”的回答是：**补上了一部分——论文已有真实电学观测的定量参照；直接支撑其独特结构机制的受控观测仍缺。** 本判断不等同于期刊编辑决定。
