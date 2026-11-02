# 初始薄椭球网格生成失败记录

首次采用OCC球体沿z直接缩放后再划分曲面网格，未获得有效网格或电学解。运行先出现反复的曲面自交/无效三角形；一次诊断运行在终端返回：

```text
Warning : ... elements remain invalid in surface 7
Surface 7 has self-intersections in its 1D mesh
Error : Identical points in triangulation: increase element size or Mesh.RandomFactor
Exception: Identical points in triangulation: increase element size or Mesh.RandomFactor
```

以上为实际终端信息的短摘录，不是完整日志文件。前一次关闭终端输出的尝试在检查无输出状态后停止。两次均未产生可用于比较的预测值。

最终实现先在伸长坐标中对原生球体和盒子划分网格，再对整个网格进行同一正行列式线性变换，得到原定薄椭球及单胞尺寸。在真实坐标中组装有限元，所有后续电学输出都来自这一实现。成功网格运行记录保存在 `mesh_initial.log` 与后续批次日志。
