# 实施计划

## 任务概述

本任务列表将 rSHUD 包现代化升级分解为可执行的编码任务。每个任务都是独立的、可测试的，并引用相关的需求。任务按照依赖关系组织，确保每一步都建立在前一步的基础上。

## 任务列表

- [x] 1. 建立基础设施和验证层
  - 更新包配置和依赖
  - 实现参数验证函数
  - 建立测试框架
  - _需求: 1, 2, 9, 10_

- [x] 1.1 更新 DESCRIPTION 文件
  - 将 terra (>= 1.7-0) 和 sf (>= 1.0-0) 移至 Imports
  - 完全移除 raster、sp、rgeos 依赖
  - 更新 R 版本要求为 >= 4.0.0
  - 移除 LinkingTo 中的 sp
  - _需求: 10_

- [x] 1.2 创建参数验证函数
  - 创建 `R/validators.R` 文件
  - 实现 `check_positive()` 函数
  - 实现 `check_file_exists()` 函数
  - 实现 `check_spatial_compatible()` 函数
  - 实现 `compatible_crs()` 和辅助函数
  - 编写完整的 roxygen2 文档
  - _需求: 4, 9, 16_

- [x] 1.3 建立测试框架
  - 创建 `tests/testthat/` 目录结构
  - 创建 `tests/testthat.R` 主测试文件
  - 创建测试辅助函数 `tests/testthat/helper-data.R`（只包含 terra/sf 测试数据）
  - 编写验证函数测试 `tests/testthat/test-validation.R`
  - _需求: 8_

- [x] 1.4 创建迁移文档
  - 创建 `inst/MIGRATION_GUIDE.md`
  - 说明不再支持 raster/sp
  - 提供 raster → terra 迁移示例
  - 提供 sp → sf 迁移示例
  - 列出常用函数替换表
  - 创建 `NEWS.md` 记录破坏性变更
  - _需求: 2, 13_

- [-] 2. 迁移核心 GIS 函数
  - 重构栅格转换函数
  - 重构投影函数
  - 确保性能提升
  - _需求: 1, 5, 12_

- [x] 2.1 重构 mesh_to_raster 函数
  - 创建 `R/gis_core.R` 文件
  - 实现 `mesh_to_raster()` 直接使用 terra::rast() 和 sf
  - 支持多种插值方法（idw, linear, nearest）
  - 使用 check_positive() 等验证函数
  - 拒绝 raster/sp 输入，提供清晰错误消息
  - 编写完整的 roxygen2 文档
  - 创建测试 `tests/testthat/test-gis-core.R`
  - _需求: 1, 5, 9, 12_

- [x] 2.2 重构 vector_to_raster 函数
  - 实现 `vector_to_raster()` 使用 terra 和 sf
  - 支持自动和手动分辨率设置
  - 支持多种聚合函数
  - 添加参数验证和错误处理
  - 编写文档和测试
  - _需求: 1, 5_

- [x] 2.3 重构坐标投影函数
  - 创建 `R/gis_projection.R` 文件
  - 重构 `crs.Albers()` 使用 sf
  - 重构 `crs.Lambert()` 使用 sf
  - 重构 `crs.long2utm()` 使用 sf
  - 实现 `project_coords()` 统一接口
  - 编写文档和测试
  - _需求: 1, 5_

- [x] 2.4 性能基准测试
  - 创建 `tests/benchmark/bench-raster.R`
  - 对比新旧实现的栅格操作性能
  - 对比新旧实现的矢量操作性能
  - 验证性能提升 >= 20%
  - 记录性能测试结果
  - _需求: 12_

- [x] 3. 迁移网格生成模块
  - 重构三角网格生成
  - 重构网格属性计算
  - 更新 S4 类定义
  - _需求: 1, 5, 16_

- [x] 3.1 重构网格生成核心函数
  - 创建 `R/mesh_generation.R` 文件
  - 重构 `generate_mesh()` (原 shud.triangle)
  - 使用 sf 处理边界和河流输入
  - 保持 RTriangle 用于三角剖分
  - 输出 sf 兼容格式
  - 编写文档
  - _需求: 1, 5_

- [x] 3.2 重构网格转换函数
  - 实现 `mesh_to_shapefile()` (原 sp.mesh2Shape)
  - 使用 sf 创建空间对象
  - 支持属性数据附加
  - 添加 CRS 支持
  - 编写文档和测试
  - _需求: 1, 5_

- [x] 3.3 重构网格属性计算
  - 创建 `R/mesh_domain.R` 文件
  - 重构 `calc_mesh_attributes()` (原 shud.att)
  - 使用 terra 提取栅格值到网格
  - 优化性能（向量化操作）
  - 编写文档和测试
  - _需求: 1, 5, 12, 16_

- [x] 3.4 更新 SHUD.MESH S4 类
  - 修改 `R/ModelClasses.R`
  - 更新 SHUD.MESH 类定义以支持 sf
  - 保持向后兼容性
  - 添加转换方法
  - 更新文档
  - _需求: 2, 5_

- [x] 3.5 编写网格模块测试
  - 创建 `tests/testthat/test-mesh.R`
  - 测试网格生成功能
  - 测试网格转换功能
  - 测试属性计算功能
  - 测试 S4 类兼容性
  - _需求: 8_

- [x] 4. 迁移河流处理模块
  - 重构河流网络构建
  - 重构河流分级和拓扑
  - 整合河流属性计算
  - _需求: 1, 5, 16, 18_

- [x] 4.1 重构河流分级函数
  - 创建 `R/river_processing.R` 文件
  - 重构 `calc_river_order()` (原 sp.RiverOrder)
  - 使用 sf 处理河流线
  - 实现 Strahler 分级算法
  - 编写文档和测试
  - _需求: 1, 5_

- [x] 4.2 重构河流路径和拓扑函数
  - 实现 `calc_river_path()` (原 sp.RiverPath)
  - 实现 `calc_river_downstream()` (原 sp.RiverDown)
  - 使用 sf 几何操作
  - 优化拓扑计算性能
  - 编写文档和测试
  - _需求: 1, 5, 12_

- [x] 4.3 整合河流属性计算函数
  - 创建 `R/river_network.R` 文件
  - 实现 `calc_river_properties()` 整合多个功能
  - 支持参数化选择计算哪些属性
  - 避免重复代码，复用现有函数
  - 编写文档和测试
  - _需求: 16, 18_

- [x] 4.4 重构河流网络构建函数
  - 实现 `build_river_network()` (原 shud.river)
  - 使用 sf 和 terra 处理输入
  - 调用 calc_river_properties 计算属性
  - 输出 sf 格式
  - 编写文档和测试
  - _需求: 1, 5, 16_

- [x] 4.5 更新 SHUD.RIVER S4 类
  - 修改 `R/ModelClasses.R`
  - 更新 SHUD.RIVER 类定义以支持 sf
  - 保持向后兼容性
  - 添加转换方法
  - 更新文档
  - _需求: 2, 5_

- [x] 4.6 编写河流模块测试
  - 创建 `tests/testthat/test-river.R`
  - 测试河流分级功能
  - 测试拓扑计算功能
  - 测试属性计算功能
  - 测试网络构建功能
  - _需求: 8_

- [x] 5. 更新 I/O 模块
  - 标准化文件读写接口
  - 确保格式兼容性
  - 支持新旧格式
  - _需求: 2, 6, 16_

- [x] 5.1 重构 SHUD 文件读取函数
  - 创建 `R/io_shud.R` 文件
  - 重命名函数：readmesh → read_mesh, readriv → read_river 等
  - 保留旧函数名作为废弃别名
  - 添加 .Deprecated 警告
  - 更新文档
  - _需求: 2, 6, 11_

- [x] 5.2 重构 SHUD 文件写入函数
  - 重命名函数：write.mesh → write_mesh, write.riv → write_river 等
  - 保留旧函数名作为废弃别名
  - 确保输出格式向后兼容
  - 更新文档
  - _需求: 2, 6, 11_

- [x] 5.3 重构时间序列 I/O 函数
  - 创建 `R/io_timeseries.R` 文件
  - 重命名 read.tsd → read_tsd, write.tsd → write_tsd
  - 保留旧函数名作为废弃别名
  - 保持 xts 格式不变
  - 更新文档
  - _需求: 2, 6, 11_

- [x] 5.4 重构模型输出读取函数
  - 创建 `R/io_output.R` 文件
  - 更新 `read_output()` 支持 terra 格式
  - 添加输出格式选项（modern/legacy）
  - 编写文档和测试
  - _需求: 1, 2_

- [x] 5.5 重构 NetCDF 处理函数
  - 创建 `R/io_netcdf.R` 文件
  - 更新 NetCDF 读取函数使用 terra
  - 保持 ncdf4 接口不变
  - 编写文档和测试
  - _需求: 1_

- [x] 5.6 编写 I/O 模块测试
  - 创建 `tests/testthat/test-io.R`
  - 测试所有读写函数
  - 测试向后兼容性
  - 测试废弃警告
  - 测试格式转换
  - _需求: 8_

- [x] 6. 更新可视化模块
  - 重构绘图函数支持新格式
  - 改进绘图质量
  - 保持接口兼容
  - _需求: 1, 2, 6_

- [x] 6.1 重构 2D 空间绘图函数
  - 创建 `R/plot_spatial.R` 文件
  - 重构 `plot_mesh_2d()` (原 map2d)
  - 支持 terra 和 sf 输入
  - 使用 terra::plot 或 ggplot2
  - 编写文档
  - _需求: 1, 2_

- [x] 6.2 重构地图对比函数
  - 实现 `compare_maps()` 支持新格式
  - 使用 ggplot2 或 terra::plot
  - 改进布局和标注
  - 编写文档
  - _需求: 1, 2_

- [x] 6.3 重构时间序列绘图函数
  - 创建 `R/plot_timeseries.R` 文件
  - 重命名 plot_tsd → plot_timeseries
  - 保留旧函数名作为别名
  - 使用 ggplot2 改进绘图
  - 编写文档
  - _需求: 2, 6_

- [x] 6.4 编写可视化模块测试
  - 创建 `tests/testthat/test-plot.R`
  - 测试绘图函数不报错
  - 测试输出格式正确
  - 使用 vdiffr 进行视觉回归测试（可选）
  - _需求: 8_

- [x] 7. 更新主接口函数
  - 重构 autoBuildModel
  - 支持新旧格式自动转换
  - 添加输出格式选项
  - _需求: 2, 15, 16_

- [x] 7.1 重构 autoBuildModel 主函数
  - 创建 `R/interface_main.R` 文件
  - 重命名为 `auto_build_model()`
  - 只接受 terra::SpatRaster 和 sf 输入
  - 拒绝 raster/sp 输入，提供迁移指导
  - 直接使用 terra/sf 函数处理
  - 返回 terra/sf 格式结果
  - 保留旧函数名作为废弃别名
  - _需求: 2, 9, 15_

- [x] 7.2 添加进度显示和日志
  - 在 auto_build_model 中添加 verbose 参数
  - 显示各阶段进度信息
  - 记录关键操作到日志
  - _需求: 15_

- [x] 7.3 实现快速建模函数
  - 实现 `quick_model()` 使用默认参数
  - 简化常见用例的接口
  - 编写文档和示例
  - _需求: 17_

- [x] 7.4 编写主接口集成测试
  - 创建 `tests/testthat/test-integration.R`
  - 测试完整的建模工作流
  - 测试 terra/sf 格式输入
  - 测试拒绝 raster/sp 输入
  - 使用小规模测试数据
  - _需求: 8, 15_

- [x] 8. 函数清理和标准化
  - 标记废弃函数
  - 移除冗余代码
  - 统一命名规范
  - _需求: 6, 11, 14, 16_

- [x] 8.1 标记所有废弃函数
  - 为所有重命名的函数添加 .Deprecated 调用
  - 在文档中标记 @deprecated
  - 在 NEWS.md 中记录所有废弃
  - _需求: 11_

- [x] 8.2 移除重复功能的函数
  - 废弃 read.df 和 write.df（推荐使用 read.table/write.table）
  - 废弃 count()（推荐使用 table()）
  - 在 NEWS.md 中记录移除
  - _需求: 14, 16_

- [x] 8.3 清理未完成或实验性函数
  - 评估 grid.subset, highlight_id, png.control 等
  - 完成或移除这些函数
  - 更新文档说明
  - _需求: 14_

- [x] 8.4 统一函数命名
  - 确保所有新函数使用 snake_case
  - 确保前缀一致（read_, write_, calc_, plot_）
  - 更新所有文档
  - _需求: 3, 6_

- [x] 9. 完善文档和测试
  - 更新所有函数文档
  - 提高测试覆盖率
  - 创建用户指南
  - _需求: 7, 8, 13_

- [x] 9.1 更新所有函数文档
  - 确保所有导出函数有完整 roxygen2 文档
  - 添加参数说明、返回值、示例
  - 添加迁移说明（@section Migration Note）
  - 添加函数分组（@family）
  - _需求: 7_

- [x] 9.2 提高测试覆盖率
  - 运行 covr::package_coverage()
  - 为覆盖率低的函数添加测试
  - 确保核心功能覆盖率 >= 70%
  - _需求: 8_

- [x] 9.3 创建迁移指南 vignette
  - 创建 `vignettes/migration-guide.Rmd`
  - 说明 v2.1 到 v2.2 的主要变更
  - 提供前后代码对比
  - 列出所有重命名的函数
  - _需求: 13_

- [x] 9.4 创建快速入门 vignette
  - 创建 `vignettes/getting-started.Rmd`
  - 展示基本工作流
  - 使用新的函数名和接口
  - 包含完整示例
  - _需求: 7_

- [x] 9.5 创建其他 vignettes
  - 创建 `vignettes/model-building.Rmd`
  - 创建 `vignettes/gis-processing.Rmd`
  - 创建 `vignettes/hydrological-analysis.Rmd`
  - _需求: 7_

- [x] 10. 最终验证和发布准备
  - 运行完整测试套件
  - 性能验证
  - 更新 NEWS 和 README
  - _需求: 所有_

- [x] 10.1 运行完整测试和检查
  - 运行 devtools::test()
  - 运行 devtools::check()
  - 确保 0 错误 0 警告
  - 修复所有问题
  - _需求: 3, 8_

- [x] 10.2 性能验证
  - 运行所有性能基准测试
  - 验证性能提升 >= 20%
  - 记录性能测试结果
  - _需求: 12_

- [x] 10.3 更新 NEWS.md
  - 记录所有主要变更
  - 列出所有废弃函数
  - 列出所有移除函数
  - 说明破坏性变更
  - 提供迁移指导
  - _需求: 11, 13_

- [x] 10.4 更新 README.md
  - 更新安装说明
  - 更新示例代码（使用新函数名）
  - 添加迁移指南链接
  - 更新徽章和链接
  - _需求: 7, 13_

- [x] 10.5 更新版本号
  - 将版本号更新为 2.2.0
  - 更新 Date 字段
  - 提交所有更改
  - 创建 git tag v2.2.0
  - _需求: 所有_

- [x] 10.6 准备发布材料
  - 构建包 tarball
  - 生成 PDF 文档
  - 准备发布说明
  - 准备 CRAN 提交材料（如适用）
  - _需求: 所有_

## 任务执行说明

### 任务标记说明
- `[ ]` - 未开始
- `[x]` - 已完成
- `*` - 可选任务（主要是额外的测试和文档）

### 执行顺序
任务按照依赖关系组织，建议按顺序执行：
1. 必须先完成父任务的所有子任务
2. 可选任务（标记 `*`）可以跳过或稍后完成
3. 每完成一个任务，运行相关测试验证

### 测试要求
- 核心实现任务必须包含单元测试
- 可选的测试任务可以在时间允许时完成
- 最终必须达到 70% 的测试覆盖率

### 文档要求
- 所有新函数和重构函数必须有完整的 roxygen2 文档
- 文档必须包含：标题、描述、参数、返回值、示例
- 迁移的函数必须包含迁移说明

### 代码审查检查点
建议在以下阶段进行代码审查：
- 完成阶段 1（基础设施）
- 完成阶段 2-4（核心功能）
- 完成阶段 7（主接口）
- 完成阶段 10（最终验证）

## 预计时间

| 阶段 | 任务 | 预计时间 |
|-----|------|---------|
| 1 | 基础设施 | 1-2 周 |
| 2 | 核心 GIS | 1-2 周 |
| 3 | 网格模块 | 1-2 周 |
| 4 | 河流模块 | 1-2 周 |
| 5 | I/O 模块 | 1 周 |
| 6 | 可视化 | 1 周 |
| 7 | 主接口 | 1 周 |
| 8 | 清理 | 1 周 |
| 9 | 文档测试 | 1 周 |
| 10 | 验证发布 | 1 周 |
| **总计** | | **9-11 周** |

## 成功标准

完成所有任务后，项目应满足：
- ✅ 所有核心函数直接使用 terra/sf（无封装层）
- ✅ 完全移除 raster/sp/rgeos 依赖
- ✅ 函数命名向后兼容（废弃别名）
- ✅ 拒绝 raster/sp 输入并提供清晰迁移指导
- ✅ 性能提升 >= 20%
- ✅ 测试覆盖率 >= 70%
- ✅ R CMD check 通过（0 错误 0 警告）
- ✅ 所有导出函数有完整文档
- ✅ 迁移指南完成（inst/MIGRATION_GUIDE.md）
- ✅ 所有 demo 使用 terra/sf 正常运行
