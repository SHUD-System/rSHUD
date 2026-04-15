# rSHUD 迁移检查清单

## 受影响的文件分析

基于代码扫描，以下文件需要进行依赖包迁移：

### 高影响文件（大量使用旧包）

#### 1. `R/Func_GIS.R` ⚠️⚠️⚠️
**使用情况**：
- `raster::extent()`, `raster::raster()`, `raster::mask()`, `raster::rasterToPoints()`
- `sp::gridded()`, `sp::fullgrid()`
- `gstat::idw()` (保留)

**迁移任务**：
- [ ] `sp2raster()` - 矢量转栅格
- [ ] `MeshData2Raster()` - 网格数据转栅格（核心函数）
- [ ] `extractRaster()` - 栅格提取
- [ ] `ForcingCoverage()` - 强迫数据覆盖

**优先级**：🔴 最高（核心功能）

#### 2. `R/GIS_RiverProcess.R` ⚠️⚠️⚠️
**使用情况**：
- `raster::extent()` - 获取空间范围
- `rgeos::gSimplify()` - 几何简化
- `sp::SpatialLines()`, `sp::Lines()`, `sp::Line()` - 空间线对象

**迁移任务**：
- [ ] `sp.RiverDown()` - 河流下游判定
- [ ] `sp.RiverPath()` - 河流路径
- [ ] `sp.RiverOrder()` - 河流分级
- [ ] `FromToNode()` - 节点连接关系

**优先级**：🔴 最高（核心功能）

#### 3. `R/River.R` ⚠️⚠️
**使用情况**：
- `sp::SpatialLinesDataFrame()` - 创建空间线数据框
- `rgdal::readOGR()` - 读取 shapefile

**迁移任务**：
- [ ] `shud.river()` - 河流网络构建
- [ ] `sp.riv2shp()` - 河流转 shapefile
- [ ] `correctRiverSlope()` - 河流坡度校正

**优先级**：🔴 最高（核心功能）

#### 4. `R/trianglulate.R` ⚠️⚠️
**使用情况**：
- `sp::SpatialPolygonsDataFrame()` - 创建空间多边形
- `rgeos::gArea()` - 计算面积
- `raster::crs()` - 坐标系统

**迁移任务**：
- [ ] `shud.triangle()` - 三角网格生成
- [ ] `sp.mesh2Shape()` - 网格转 shapefile

**优先级**：🔴 最高（核心功能）

#### 5. `R/MeshDomain.R` ⚠️⚠️
**使用情况**：
- `apply.raster()` - 栅格数据应用到网格

**迁移任务**：
- [ ] `shud.att()` - 网格属性生成
- [ ] `apply.raster()` - 辅助函数

**优先级**：🔴 最高（核心功能）

#### 6. `R/GIS_delineation.R` ⚠️⚠️
**使用情况**：
- `raster::raster()` - 读取栅格
- `raster::extent()`, `raster::crs()` - 栅格属性
- `rgdal::readOGR()` - 读取矢量

**迁移任务**：
- [ ] `watershedDelineation()` - 流域划分

**优先级**：🟡 中等（依赖 whitebox，需要测试）

### 中影响文件

#### 7. `R/plotMap.R` ⚠️
**使用情况**：
- `raster::plot()`, `raster::contour()`, `raster::stack()`
- `is.Raster()` - 类型检查

**迁移任务**：
- [ ] `compareMaps()` - 多图对比
- [ ] `plot_animate()` - 动画绘制
- [ ] `plot_sp()` - 空间数据绘图

**优先级**：🟡 中等（可视化功能）

#### 8. `R/plot.R` ⚠️
**使用情况**：
- `raster::plot()` - 栅格绘图
- `MeshData2Raster()` - 调用核心函数

**迁移任务**：
- [ ] `map2d()` - 2D 地图绘制

**优先级**：🟡 中等（可视化功能）

#### 9. `R/readout.R` ⚠️
**使用情况**：
- `raster::focal()` - 焦点统计
- `MeshData2Raster()` - 调用核心函数

**迁移任务**：
- [ ] `readout()` - 读取输出文件

**优先级**：🟡 中等（I/O 功能）

#### 10. `R/NetCDF.R` ⚠️
**使用情况**：
- `raster::raster()`, `raster::extent()`, `raster::stack()`
- `xyz2Raster()` - 坐标转栅格

**迁移任务**：
- [ ] `readnc()` - 读取 NetCDF
- [ ] `xyz2Raster()` - 辅助函数

**优先级**：🟡 中等（I/O 功能）

#### 11. `R/GIS_Projection.R` ⚠️
**使用情况**：
- `sp::CRS()`, `sp::spTransform()` - 坐标转换
- `raster::extent()` - 范围提取

**迁移任务**：
- [ ] `crs.Albers()` - Albers 投影
- [ ] `crs.Lambert()` - Lambert 投影
- [ ] `crs.long2utm()` - 经纬度转 UTM

**优先级**：🟡 中等（GIS 工具）

#### 12. `R/GIS_SimplifybyLen.R` ⚠️
**使用情况**：
- `methods::as(sp, 'SpatialLines')` - 类型转换

**迁移任务**：
- [ ] `sp.simplifyLen()` - 按长度简化
- [ ] `SimplifybyLen()` - 简化函数

**优先级**：🟢 低（辅助功能）

#### 13. `R/autoBuildModel.R` ⚠️
**使用情况**：
- `is(spr, 'SpatialLines')` - 类型检查
- `sp::SpatialLinesDataFrame()` - 创建空间数据

**迁移任务**：
- [ ] `autoBuildModel()` - 自动建模主函数

**优先级**：🔴 最高（主要接口）

### 不受影响的文件 ✓

以下文件不直接使用 raster/sp/rgeos：
- `R/DataFilter.R` - 数据过滤
- `R/Func_landcover.R` - 土地覆盖
- `R/Func_Misc.R` - 杂项函数
- `R/Func_PTF.R` - 土壤转换函数
- `R/Hydro_obs.R` - 水文观测
- `R/ModelBC.R` - 边界条件
- `R/ModelClasses.R` - 类定义
- `R/ModelConfigure.R` - 模型配置
- `R/ModelInfo.R` - 模型信息
- `R/NLCD.R` - NLCD 数据
- `R/PET.R` - 蒸散发
- `R/Project.R` - 项目管理
- `R/RcppExports.R` - Rcpp 导出
- `R/readinput.R` - 输入读取（文本格式）
- `R/readTSD.R` - 时间序列读取
- `R/SHUD_Env.R` - 环境设置
- `R/Sinks.R` - 汇点处理
- `R/WaterBalance.R` - 水量平衡
- `R/writeInput.R` - 输入写入

## 迁移优先级总结

### 第一批（核心 GIS 和网格功能）
1. `Func_GIS.R` - 特别是 `MeshData2Raster()`
2. `trianglulate.R` - 网格生成
3. `MeshDomain.R` - 网格属性
4. `GIS_RiverProcess.R` - 河流处理
5. `River.R` - 河流网络

### 第二批（主接口和投影）
6. `autoBuildModel.R` - 主函数
7. `GIS_Projection.R` - 坐标投影
8. `GIS_delineation.R` - 流域划分

### 第三批（I/O 和可视化）
9. `NetCDF.R` - NetCDF 处理
10. `readout.R` - 输出读取
11. `plot.R` - 基础绘图
12. `plotMap.R` - 地图绘制

### 第四批（辅助功能）
13. `GIS_SimplifybyLen.R` - 简化工具

## 关键函数依赖关系

```
autoBuildModel()
├── shud.triangle()
│   ├── sp.mesh2Shape() [需要迁移]
│   └── MeshData2Raster() [需要迁移]
├── shud.river()
│   ├── sp.RiverOrder() [需要迁移]
│   ├── sp.RiverPath() [需要迁移]
│   └── sp.riv2shp() [需要迁移]
├── shud.att()
│   └── apply.raster() [需要迁移]
└── watershedDelineation()
    └── [依赖 whitebox + raster]
```

## 迁移策略建议

### 方案 A：渐进式迁移（推荐）
1. 创建新的辅助函数支持 terra/sf
2. 保留旧函数，添加 `.Deprecated()` 警告
3. 逐步迁移核心函数
4. 更新文档和示例
5. 发布 v3.0.0

**优点**：风险低，用户有时间适应
**缺点**：维护成本高，代码冗余

### 方案 B：一次性迁移
1. 直接替换所有旧包调用
2. 更新所有依赖函数
3. 发布 v3.0.0（破坏性更新）

**优点**：代码简洁，维护成本低
**缺点**：风险高，可能破坏现有工作流

### 推荐：混合方案
- 核心函数（第一批）：完全迁移到 terra/sf
- 提供转换函数：`as_terra()`, `as_sf()` 用于兼容旧对象
- 主接口函数：自动检测输入类型，支持新旧格式
- 文档：明确标注新旧 API 差异

## 下一步

1. ✅ 完成迁移清单
2. ⏳ 创建兼容层函数
3. ⏳ 开始第一批核心函数迁移
4. ⏳ 编写单元测试
5. ⏳ 更新文档
