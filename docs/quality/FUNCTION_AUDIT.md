# rSHUD 函数审查和清理计划

## 审查标准

每个函数按以下标准评估：

- **使用频率**：高/中/低
- **功能完整性**：完整/部分/未完成
- **文档质量**：完整/部分/缺失
- **测试覆盖**：有/无
- **代码质量**：好/一般/差
- **决策**：保留/重构/废弃/移除

## 函数分类

### A 类：核心函数（必须保留和优化）

#### 模型构建
| 函数 | 文件 | 状态 | 优先级 | 行动 |
|------|------|------|--------|------|
| `autoBuildModel()` | autoBuildModel.R | 核心 | 🔴 最高 | 重构，迁移到 terra/sf |
| `shud.triangle()` | trianglulate.R | 核心 | 🔴 最高 | 重构，改进文档 |
| `shud.mesh()` | MeshDomain.R | 核心 | 🔴 最高 | 重构 |
| `shud.river()` | River.R | 核心 | 🔴 最高 | 重构 |
| `shud.att()` | MeshDomain.R | 核心 | 🔴 最高 | 重构 |

#### I/O 函数
| 函数 | 文件 | 状态 | 优先级 | 行动 |
|------|------|------|--------|------|
| `readmesh()` | readinput.R | 完整 | 🟡 中 | 重命名为 `read_mesh()` |
| `readriv()` | readinput.R | 完整 | 🟡 中 | 重命名为 `read_river()` |
| `readatt()` | readinput.R | 完整 | 🟡 中 | 重命名为 `read_att()` |
| `readpara()` | readinput.R | 完整 | 🟡 中 | 重命名为 `read_para()` |
| `readcalib()` | readinput.R | 完整 | 🟡 中 | 重命名为 `read_calib()` |
| `readic()` | readinput.R | 完整 | 🟡 中 | 重命名为 `read_ic()` |
| `write.mesh()` | writeInput.R | 完整 | 🟡 中 | 重命名为 `write_mesh()` |
| `write.riv()` | writeInput.R | 完整 | 🟡 中 | 重命名为 `write_river()` |
| `readout()` | readout.R | 完整 | 🟡 中 | 重构，支持 terra |
| `read.tsd()` | readTSD.R | 完整 | 🟡 中 | 重命名为 `read_tsd()` |
| `write.tsd()` | readTSD.R | 完整 | 🟡 中 | 重命名为 `write_tsd()` |

#### GIS 处理
| 函数 | 文件 | 状态 | 优先级 | 行动 |
|------|------|------|--------|------|
| `MeshData2Raster()` | Func_GIS.R | 核心 | 🔴 最高 | 重构为 `mesh_to_raster()` |
| `sp2raster()` | Func_GIS.R | 常用 | 🔴 最高 | 重构为 `vector_to_raster()` |
| `watershedDelineation()` | GIS_delineation.R | 重要 | 🟡 中 | 重构为 `delineate_watershed()` |
| `sp.RiverOrder()` | GIS_RiverProcess.R | 核心 | 🔴 最高 | 重构为 `river_order()` |
| `sp.RiverPath()` | GIS_RiverProcess.R | 核心 | 🔴 最高 | 重构为 `river_path()` |
| `sp.RiverDown()` | GIS_RiverProcess.R | 核心 | 🔴 最高 | 重构为 `river_downstream()` |

### B 类：重要函数（需要改进）

#### 水文计算
| 函数 | 文件 | 状态 | 优先级 | 行动 |
|------|------|------|--------|------|
| `PET_PM()` | PET.R | 完整 | 🟡 中 | 改进文档，添加测试 |
| `MeltFactor()` | PET.R | 完整 | 🟡 中 | 改进文档 |
| `wb.all()` | WaterBalance.R | 完整 | 🟡 中 | 重命名为 `water_balance()` |
| `wb.ele()` | WaterBalance.R | 完整 | 🟡 中 | 重命名为 `wb_element()` |
| `wb.riv()` | WaterBalance.R | 完整 | 🟡 中 | 重命名为 `wb_river()` |

#### 可视化
| 函数 | 文件 | 状态 | 优先级 | 行动 |
|------|------|------|--------|------|
| `map2d()` | plot.R | 常用 | 🟡 中 | 重构，支持 terra |
| `plot_tsd()` | plotMap.R | 常用 | 🟡 中 | 改进，使用 ggplot2 |
| `hydrograph()` | Hydro_obs.R | 常用 | 🟡 中 | 改进 |
| `compareMaps()` | plotMap.R | 有用 | 🟢 低 | 重构 |

#### 数据处理
| 函数 | 文件 | 状态 | 优先级 | 行动 |
|------|------|------|--------|------|
| `extractCoords()` | ModelConfigure.R | 常用 | 🟡 中 | 移到 compat 层 |
| `ProjectCoordinate()` | GIS_Projection.R | 常用 | 🟡 中 | 重构为 `project_coords()` |
| `crs.Albers()` | GIS_Projection.R | 有用 | 🟡 中 | 保留，改进文档 |
| `crs.Lambert()` | GIS_Projection.R | 有用 | 🟡 中 | 保留，改进文档 |

### C 类：辅助函数（评估后决定）

#### 保留但改进
| 函数 | 文件 | 问题 | 行动 |
|------|------|------|------|
| `SimplifybyLen()` | GIS_SimplifybyLen.R | 命名不规范 | 重命名为 `simplify_by_length()` |
| `correctRiverSlope()` | River.R | 文档不足 | 改进文档和错误处理 |
| `RiverLength()` | River.R | 功能简单 | 考虑合并到 `river_properties()` |
| `RiverSlope()` | River.R | 功能简单 | 考虑合并到 `river_properties()` |
| `fishnet()` | Func_GIS.R | 有用 | 改进文档 |
| `voronoipolygons()` | Func_GIS.R | 有用 | 改进文档 |

#### 需要评估
| 函数 | 文件 | 问题 | 建议 |
|------|------|------|------|
| `plot_animate()` | plotMap.R | 依赖复杂，使用少 | 考虑移到单独包或示例 |
| `ts2map()` | plotMap.R | 功能重复 | 考虑合并到其他函数 |
| `meshSinks()` | Sinks.R | 文档缺失 | 补充文档或移除 |
| `datafilter.riv()` | DataFilter.R | 功能不清晰 | 重新设计或移除 |

### D 类：考虑废弃或移除

#### 功能重复
| 函数 | 文件 | 原因 | 行动 |
|------|------|------|------|
| `read.df()` | readinput.R | 与 `read.table()` 功能重复 | 废弃，推荐使用标准函数 |
| `write.df()` | writeInput.R | 与 `write.table()` 功能重复 | 废弃，推荐使用标准函数 |
| `count()` | Func_Misc.R | 与 `table()` 功能重复 | 移除 |

#### 未完成或实验性
| 函数 | 文件 | 原因 | 行动 |
|------|------|------|------|
| `grid.subset()` | Func_GIS.R | 功能不完整 | 完成或移除 |
| `highlight_id()` | plotMap.R | 使用场景有限 | 移到示例代码 |
| `png.control()` | plotMap.R | 功能过于简单 | 移除，用户可直接使用 `png()` |

#### 过时或不再需要
| 函数 | 文件 | 原因 | 行动 |
|------|------|------|------|
| `BasicPlot()` | plot.R | 被更好的函数替代 | 废弃 |
| `SimpleSpatial()` | Func_GIS.R | 功能不明确 | 评估后决定 |

## 清理计划

### 阶段 1：标记废弃函数（v3.0）
```r
# 添加废弃警告
read.df <- function(...) {
  .Deprecated("read.table", 
              msg = "read.df() is deprecated. Use read.table() or readr::read_delim() instead.")
  read.table(...)
}
```

### 阶段 2：重命名核心函数（v3.0）
```r
# 新函数
read_mesh <- function(...) {
  # 实现
}

# 保留旧函数作为别名
readmesh <- function(...) {
  .Deprecated("read_mesh")
  read_mesh(...)
}
```

### 阶段 3：移除废弃函数（v3.1 或 v4.0）
- 在 v3.0 中标记为废弃
- 在 v3.1 中移除（或等到 v4.0）
- 在 NEWS.md 中明确说明

## 新增函数建议

### 高级接口函数
```r
# 统一的河流属性计算
river_properties <- function(river, dem, ...) {
  list(
    order = river_order(river),
    length = river_length(river),
    slope = river_slope(river, dem),
    downstream = river_downstream(river)
  )
}

# 统一的网格属性计算
mesh_properties <- function(mesh, ...) {
  list(
    area = mesh_area(mesh),
    centroid = mesh_centroid(mesh),
    neighbors = mesh_neighbors(mesh)
  )
}
```

### 便捷函数
```r
# 快速模型设置
quick_model <- function(domain, dem, output_dir, ...) {
  # 使用默认参数快速建立模型
}

# 批量处理
batch_process <- function(models, fun, ...) {
  # 批量处理多个模型
}
```

### 诊断函数
```r
# 模型诊断
diagnose_model <- function(model) {
  # 检查模型配置的常见问题
}

# 数据质量检查
check_data_quality <- function(data, type = c("mesh", "river", "forcing")) {
  # 检查输入数据质量
}
```

## 文档改进计划

### 创建主题文档（Vignettes）
1. **Getting Started** - 快速入门
2. **Model Building** - 模型构建详解
3. **GIS Processing** - GIS 数据处理
4. **Hydrological Analysis** - 水文分析
5. **Migration Guide** - v2.x 到 v3.0 迁移指南

### 改进函数分组
```r
# 在文档中按功能分组
#' @family model-building
#' @family io-functions
#' @family gis-processing
#' @family hydrological-analysis
#' @family visualization
```

## 实施时间表

### 第 1-2 周：审查和规划
- [x] 完成函数审查
- [ ] 确定优先级
- [ ] 制定详细计划

### 第 3-4 周：核心函数重构
- [ ] 迁移 A 类函数到 terra/sf
- [ ] 添加参数验证
- [ ] 改进错误处理

### 第 5-6 周：重命名和标准化
- [ ] 重命名函数（保留别名）
- [ ] 统一命名规范
- [ ] 更新文档

### 第 7 周：清理和优化
- [ ] 标记废弃函数
- [ ] 移除冗余代码
- [ ] 性能优化

### 第 8 周：测试和文档
- [ ] 添加单元测试
- [ ] 更新所有文档
- [ ] 创建 vignettes

## 成功指标

- [ ] 所有 A 类函数迁移到 terra/sf
- [ ] 所有导出函数有完整文档
- [ ] 测试覆盖率 > 70%
- [ ] 通过 `R CMD check` 无错误无警告
- [ ] 所有 demo 脚本正常运行
- [ ] 性能提升 > 20%（对于大数据集）

## 风险和缓解

### 风险 1：破坏现有工作流
**缓解**：
- 保留旧函数作为别名
- 提供详细的迁移指南
- 在多个版本中逐步废弃

### 风险 2：测试不充分
**缓解**：
- 使用现有 demo 作为集成测试
- 添加单元测试
- Beta 测试阶段征求用户反馈

### 风险 3：文档更新不及时
**缓解**：
- 在重构时同步更新文档
- 使用 roxygen2 自动生成
- 代码审查时检查文档

### 风险 4：性能退化
**缓解**：
- 性能基准测试
- 优化关键路径
- 使用 profiling 工具识别瓶颈
