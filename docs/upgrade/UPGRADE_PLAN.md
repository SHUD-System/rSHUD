# rSHUD 包升级计划

## 升级目标

1. **依赖包现代化**：从过时的 `raster`/`sp`/`rgeos` 迁移到 `terra`/`sf`
2. **提高容错性和灵活性**：增强函数的参数验证和错误处理
3. **代码标准化**：符合 R 包开发规范（tidyverse style guide）
4. **功能重构**：重新评估和优化现有函数
5. **清理冗余代码**：移除未完成或用处有限的函数

## 当前依赖分析

### 需要替换的包
- `raster` → `terra` (现代化的栅格处理)
- `sp` → `sf` (现代化的矢量处理)
- `rgeos` → `sf` (几何操作已集成在 sf 中)
- `rgdal` → `sf` (I/O 已集成在 sf 中)

### 保留的包
- `Rcpp` - C++ 集成
- `ggplot2` - 可视化
- `xts`, `zoo` - 时间序列
- `hydroGOF` - 水文指标
- `gstat` - 地统计
- `RTriangle` - 三角网格生成
- `ncdf4` - NetCDF 处理

## 升级策略

### 阶段 1：依赖包迁移准备（1-2周）
- [ ] 创建兼容层函数，支持新旧两种格式
- [ ] 更新 DESCRIPTION 文件
- [ ] 创建迁移指南文档

### 阶段 2：核心函数重构（2-3周）
- [ ] GIS 处理函数（`Func_GIS.R`, `GIS_*.R`）
- [ ] 网格生成函数（`trianglulate.R`, `MeshDomain.R`）
- [ ] 河流处理函数（`River.R`, `GIS_RiverProcess.R`）
- [ ] I/O 函数（`readinput.R`, `writeInput.R`）

### 阶段 3：可视化和分析函数（1-2周）
- [ ] 绘图函数（`plot.R`, `plotMap.R`）
- [ ] 水文分析函数（`WaterBalance.R`, `Hydro_obs.R`）

### 阶段 4：清理和优化（1周）
- [ ] 移除冗余函数
- [ ] 统一命名规范
- [ ] 增强错误处理
- [ ] 添加单元测试

### 阶段 5：文档和示例更新（1周）
- [ ] 更新所有函数文档
- [ ] 更新 demo 脚本
- [ ] 创建迁移指南
- [ ] 更新 README

## 关键迁移映射

### raster → terra
```r
# 旧代码
r <- raster::raster(file)
s <- raster::stack(files)
e <- raster::extent(r)
raster::projection(r) <- crs

# 新代码
r <- terra::rast(file)
s <- terra::rast(files)  # stack 和 brick 统一为 rast
e <- terra::ext(r)
terra::crs(r) <- crs
```

### sp → sf
```r
# 旧代码
sp <- sp::SpatialPolygonsDataFrame(...)
sp <- rgeos::gSimplify(sp, tol)
area <- rgeos::gArea(sp)

# 新代码
sf <- sf::st_as_sf(...)
sf <- sf::st_simplify(sf, dTolerance = tol)
area <- sf::st_area(sf)
```

## 函数审查清单

### 高优先级（核心功能）
- ✓ `autoBuildModel()` - 自动建模主函数
- ✓ `shud.triangle()` - 网格生成
- ✓ `shud.river()` - 河流网络
- ✓ `read*()` / `write*()` - I/O 函数
- ✓ `MeshData2Raster()` - 网格数据转换

### 中优先级（常用功能）
- ○ `watershedDelineation()` - 流域划分
- ○ `PET_PM()` - 蒸散发计算
- ○ `map2d()` - 2D 可视化
- ○ `wb.*()` - 水量平衡

### 低优先级（辅助功能）
- △ `plot_animate()` - 动画绘制
- △ 部分 GIS 辅助函数

### 待评估/可能移除
- ✗ 未完成的实验性函数
- ✗ 重复功能的函数
- ✗ 依赖已废弃包且无替代的函数

## 代码规范要求

1. **命名规范**
   - 函数名：snake_case 或 camelCase（保持一致）
   - 变量名：snake_case
   - 常量：UPPER_CASE

2. **参数验证**
   ```r
   function(x, y = NULL) {
     # 参数检查
     if (!inherits(x, "SpatRaster")) {
       stop("x must be a SpatRaster object")
     }
     if (is.null(y)) {
       y <- default_value()
     }
     # 函数逻辑
   }
   ```

3. **错误处理**
   - 使用 `stop()` 报告错误
   - 使用 `warning()` 报告警告
   - 使用 `message()` 报告信息

4. **文档规范**
   - 所有导出函数必须有完整的 roxygen2 文档
   - 包含示例代码
   - 说明参数类型和返回值

## 测试策略

1. 为核心函数添加单元测试
2. 使用示例数据集进行集成测试
3. 确保向后兼容性（至少一个版本）

## 风险控制

1. **版本管理**
   - 主版本号升级：3.0.0（破坏性变更）
   - 提供迁移指南

2. **兼容性**
   - 保留旧函数作为 deprecated（标记为废弃）
   - 提供 6-12 个月的过渡期

3. **回归测试**
   - 使用现有 demo 作为测试用例
   - 确保输出结果一致

## 下一步行动

1. 审查并确认升级计划
2. 创建开发分支 `dev-v3.0`
3. 开始阶段 1：依赖包迁移准备
