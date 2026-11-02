# 公开岩石电学观测的外部数值比较

本目录对Han等（2022）人工砂岩的公开电学观测执行固定输入的三维有限元预测。研究状态固定在有效压力5 MPa；I-0/I-90校准背景，I-45独立检验完整岩心方向响应，F-0/F-45/F-90为含裂缝预测目标。六根岩心的其余压力记录保存在数据表中，没有增加压力扫描。

先阅读 `REPORT_CN.md` 中的实际结果和 `VALIDATION_ASSESSMENT.md` 中的证据判断。`DATA_SOURCES.md`列出测量、几何、物性、假设及取得方式。

## 文件与目录

| 文件/目录 | 内容 |
|---|---|
| `REPORT_CN.md` | 中文结果报告、逐项目标比较和数值依据 |
| `VALIDATION_ASSESSMENT.md` | 外部比较支持的内容、尚未直接实测的主命题 |
| `DATA_SOURCES.md` | DOI、表格、参数出处及实际数据取得情况 |
| `INPUT_PROTOCOL.json` | 首轮求解前记录的校准/目标分离与名义模型；增加第四档网格见下项 |
| `RUN_AMENDMENTS.md` | 网格生成修复、第四档加密和容差传播的实际执行记录 |
| `EXPERIMENT_DATA.csv` | 六根岩心×十档压力的60行原始观测转录 |
| `MODEL_PREDICTIONS.csv` | 每个已算场景、样品和网格的绝对预测电导率及背景值 |
| `COMPARISON.csv` | 逐项实测值、模型值、带符号误差、绝对误差及增量估计 |
| `MEASUREMENT_PROPAGATION.csv` | 原文电导率/孔隙率容差的端点传播结果 |
| `UNCERTAINTY_BUDGET.csv` | 固定名义模型下的观测和预测容差区间 |
| `VALIDATION_FIGURE.pdf` | 从真实预测CSV生成的实验—模型对照图 |
| `scripts/` | 贴体网格、有限元、主批、敏感性、容差传播及交付生成脚本 |
| `runs/` | 网格、单胞/圆柱求解缓存、名义与敏感性明细、I-45独立检验 |
| `checks/` | 有限元解析对照及输出核验记录 |
| `sources/` | 参数来源说明、独立数值转录及原核对结果；原文与手册通过官方链接访问 |
| `data/` | 已与独立转录核对的公开Han2022表1数值输入；不包含用户任务方案或通信草稿 |
| `runs_*.log`, `mesh_initial.log` | 原始运行输出和网格生成记录 |

## 计算环境

本次实际运行环境为Windows、Python3.10.20：

| 包 | 版本 |
|---|---:|
| gmsh | 4.15.2 |
| numpy | 2.2.6 |
| scipy | 1.15.3 |
| pandas | 2.3.3 |
| matplotlib | 3.10.9 |
| pyamg | 5.4.1 |

本机解释器为 `D:/Users/PC/anaconda3/envs/efield310/python.exe`。复算可使用同版本的独立Python环境，无需COMSOL或MATLAB。字体绘制优先使用Times New Roman；缺少该字体不改变CSV和数值结果。

建立环境的示例命令：

```powershell
py -3.10 -m venv .venv
.\.venv\Scripts\python.exe -m pip install gmsh==4.15.2 numpy==2.2.6 scipy==1.15.3 pandas==2.3.3 matplotlib==3.10.9 pyamg==5.4.1
```

以下命令中的`python`应指向上述环境。先将工作目录切换到本目录。

## 复算顺序

### 1. 有限元解析检查

```powershell
python scripts/tet_fem.py checks/FEM_NUMERICAL_CHECKS.json
```

检查包含均匀标量/张量导电、倾斜四面体梯度与能量、两相串联/并联界面传输。完整端面电流由Dirichlet反力求和，体积分功率用于独立核对。

### 2. 名义预测及四档网格

```powershell
python scripts/run_validation.py --stage nominal
```

主批为F-0/F-45/F-90×四档网格level0–3。另计算I-45在30/40/50 mm长度下的四档网格检验，以及±45°对称性。生成 `runs/predictions_nominal.csv`、`runs/intact_holdout.csv`、`runs/operator_symmetry_check.json` 并更新总预测表及比较表。

### 3. 固定情景敏感性

```powershell
python scripts/run_validation.py --stage sensitivity
```

逐项改变开度比例0/0.5、层距2.5/3.5 mm、盐水电导率5.2 S/m、背景孔隙率归一化及圆柱长度30/50 mm。其余参数保留名义值；情景采用level1网格，与该档名义解对照。基质域孔隙率闭合另外计算level3，供两种闭合在同一最细网格比较。生成 `runs/predictions_sensitivity.csv` 并更新总表。该流程没有按实验误差选取“最佳”场景。

### 4. 测量容差传播

```powershell
python scripts/uncertainty_budget.py
```

使用原文电导率±2%、孔隙率±0.0001，将I-0/I-90背景校准误差和F孔隙率误差传播到固定名义几何的预测。使用level2网格。输出是确定性容差包络，与几何/物性情景分开记录。

### 5. 生成图、数字摘要与最终核验

```powershell
python scripts/build_deliverables.py
python scripts/validate_outputs.py
```

交付生成程序读取已完成求解的CSV，生成图件、数字摘要和环境记录；它不拟合或修改实验测量。中文报告根据这些实际结果撰写；论文插入草稿不在公开包内。最后的核验检查数据/预测对应、误差算式、网格与输出文件。具体核验内容和结果以 `checks/` 中实际生成的记录为准。

## 缓存与重新计算

`runs/meshes/`保存单胞和圆柱网格；`runs/cell_results/`及`runs/specimen_results/`保存已经实际求解的JSON。重复运行优先使用这些缓存。若要从头复算，应在新的工作副本中保留 `scripts/`、`data/`、`INPUT_PROTOCOL.json`等输入，使用空的`runs/`目录，再按以上顺序运行。保留交付目录中的原始运行结果，便于逐项对照。

主入口中的参数与`INPUT_PROTOCOL.json`对应；该JSON是输入记录，脚本没有把它当作通用配置解析器。改变模型参数或数值实现时，应同时更新输入记录，并在新的运行目录中求解，以免混用旧缓存。

程序顺序执行不同样品，复用已有网格并限制峰值内存；底层线性代数线程默认设为2。`.npz`网格及JSON缓存可直接检查，网格坐标存为mm，有限元组装转换为m，电导率单位为S/m。

## 模型和数据含义

实际实验是1 kHz交流观测。模型求解稳态导电方程，以该频率下报告的导电响应为比较对象。完整背景按原文方向Archie关系校准，校准只用I-0/I-90。I-45和三个F样品的电导率均未用于模型参数调整。

裂缝为有限扁椭球，位于镜面对称规则单胞中心；单胞体积按制造尺度与裂缝密度确定。为获得薄椭球贴体网格，程序先在伸长坐标中生成球—盒网格，再将节点映射到真实椭球坐标，并在真实坐标中组装P1四面体有限元。网格没有通过元素中心染色定义裂缝。

单胞主轴电导率随后组成各向异性张量，在真实有限圆柱内旋转并求完整端面等势/侧壁绝缘问题。45°结果来自此测量算子。数值圆柱横截面积采用网格体积除以长度校正多边形近似，几何面积因子单独保存。

4.69 S/m为Han2018同名盐度向Han2022的转移值；5.2 S/m为独立物性情景。开度比例、间距和长度情景是固定假设的敏感性计算。孔隙率闭合分别采用原文总体积归一化和剩余基质域归一化，并保存两种结果。

`predicted_fracture_increment_S_per_m`是同一模型含裂缝值减其背景值。`measured_minus_transferred_background_S_per_m`是实测F电导率减模型转移背景，两者的证据性质不同。实验没有提供每根F岩心加入裂缝前的同体背景测量。

本次最终判断以总量、方向、增量及数值敏感性分别给出。它提供公开真实观测下的部分外部检验；本文严格A/B背景匹配、外围结构移动、排序反转及模态预测仍需对应的受控电学实验。

## 原始来源与公开归档

Han2022表1已独立复核。Han2018的原始表S1未取得；Mohammed2020两份原始XLSX已取得但未运行对应模型，具体见 `DATA_SOURCES.md`。本目录不重打包未取得的论文全文。若公开发布本目录，请保留数据DOI、原数据许可及处理说明；第三方原文/手册按各来源的许可提供链接。


## 公开派生版本说明

本目录为v1.2.0公开派生版，保留原始数值输入、实际FE结果、网格缓存、固定敏感性和容差传播。仅更换输入表路径并去掉内部论文草稿的交付检查要求；FE数值算法和结果没有改动。`VALIDATION_FIGURE.{pdf,png,svg}`是当前Figure S1的同源图件，生成数据和入口均保留。

原私人任务方案、邮件草稿、论文插入建议、第三方全文PDF/长篇摘录及网页代码未打包；来源改用`DATA_SOURCES.md`中的DOI和官方链接。原始提取程序依赖出版商摘录，因此不随本版分发；可用独立转录CSV与输入表复核数字，重新对原表核对需要访问出版商DOI。`checks/independent_arithmetic_audit.json`是36行阶段的历史检查，最终39行及173项PDE诊断见`checks/DELIVERY_VALIDATION.json`，两者没有混称同一阶段。

代码遵循项目GPL-3.0-only；项目生成结果与图数据为CC BY 4.0。第三方测量保留原作者和DOI归属，本文档不授予第三方原文的再分发权。更多范围见`PUBLIC_ARCHIVE_NOTE.md`。
