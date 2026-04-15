# 设计文档

## 概述

本设计文档描述了 rSHUD 包现代化升级的技术方案。升级将包从传统空间库（raster/sp/rgeos）迁移到现代空间库（terra/sf），同时提高代码质量、改善文档和测试覆盖率。设计遵循最小工程化原则和渐进式实现策略。

**版本规划**：当前版本 2.1.0 → 新版本 2.2.0（次要版本更新）

### 设计目标

1. **性能提升**：通过直接使用 terra/sf 实现 20-50% 的性能改进
2. **现代化**：完全迁移到 terra/sf，移除所有传统空间库依赖
3. **简洁性**：直接使用 terra/sf API，不提供不必要的封装层
4. **代码质量**：标准化命名、完整文档、全面测试
5. **可维护性**：模块化设计、单一职责、避免重复
6. **清晰迁移**：提供详细迁移指南，清晰的错误消息

## 架构

### 整体架构

```
rSHUD v2.1.0 架构
├── 用户接口层
│   ├── 主接口函数（autoBuildModel 等）
│   └── 专用函数（read_*, write_*, calc_* 等）
├── 验证层
│   └── 参数验证函数（check_positive, check_file_exists 等）
├── 核心功能层
│   ├── GIS 处理模块
│   ├── 网格生成模块
│   ├── 河流处理模块
│   ├── I/O 模块
│   └── 水文计算模块
└── 底层库
    ├── 现代空间库（terra, sf）- 直接使用
    └── 其他依赖（Rcpp, ggplot2 等）
```

### 模块划分

#### 1. 验证模块（Validation Layer）
**文件**：`R/validators.R`

**职责**：
- 参数验证和类型检查
- 空间对象兼容性验证
- 提供清晰的错误消息

**关键函数**：
- `check_positive()` - 验证正数参数
- `check_file_exists()` - 验证文件存在
- `check_spatial_compatible()` - 验证空间对象兼容性
- `compatible_crs()` - CRS 兼容性检查
- `crs_to_string()` - CRS 转字符串（内部）
- `extract_epsg()` - 提取 EPSG 代码（内部）

#### 2. GIS 处理模块
**文件**：
- `R/gis_core.R`（新建）- 核心 GIS 函数
- `R/gis_projection.R`（重构自 GIS_Projection.R）
- `R/gis_delineation.R`（重构自 GIS_delineation.R）

**职责**：
- 栅格-矢量转换
- 坐标投影
- 流域划分
- 空间分析

**关键函数**：
- `mesh_to_raster()` - 网格数据转栅格（重构 MeshData2Raster）
- `vector_to_raster()` - 矢量转栅格（重构 sp2raster）
- `project_coords()` - 坐标投影（重构 ProjectCoordinate）
- `delineate_watershed()` - 流域划分（重构 watershedDelineation）

#### 3. 网格生成模块
**文件**：
- `R/mesh_generation.R`（重构自 trianglulate.R）
- `R/mesh_domain.R`（重构自 MeshDomain.R）
- `R/mesh_utils.R`（新建）

**职责**：
- 三角网格生成
- 网格属性计算
- 网格拓扑关系

**关键函数**：
- `generate_mesh()` - 生成三角网格（重构 shud.triangle）
- `mesh_to_shapefile()` - 网格转 shapefile（重构 sp.mesh2Shape）
- `calc_mesh_attributes()` - 计算网格属性（重构 shud.att）
- `mesh_topology()` - 网格拓扑（使用 Rcpp）

#### 4. 河流处理模块
**文件**：
- `R/river_network.R`（重构自 River.R）
- `R/river_processing.R`（重构自 GIS_RiverProcess.R）
- `R/river_utils.R`（新建）

**职责**：
- 河流网络构建
- 河流分级
- 河流拓扑关系
- 河流属性计算

**关键函数**：
- `build_river_network()` - 构建河流网络（重构 shud.river）
- `calc_river_order()` - 河流分级（重构 sp.RiverOrder）
- `calc_river_path()` - 河流路径（重构 sp.RiverPath）
- `calc_river_downstream()` - 下游关系（重构 sp.RiverDown）
- `calc_river_properties()` - 河流属性（整合多个函数）

#### 5. I/O 模块
**文件**：
- `R/io_shud.R`（重构自 readinput.R, writeInput.R）
- `R/io_timeseries.R`（重构自 readTSD.R）
- `R/io_output.R`（重构自 readout.R）
- `R/io_netcdf.R`（重构自 NetCDF.R）

**职责**：
- SHUD 格式文件读写
- 时间序列数据处理
- 模型输出读取
- NetCDF 文件处理

**关键函数**：
- `read_mesh()`, `write_mesh()` - 网格文件
- `read_river()`, `write_river()` - 河流文件
- `read_att()`, `write_att()` - 属性文件
- `read_tsd()`, `write_tsd()` - 时间序列
- `read_output()` - 模型输出

#### 6. 水文计算模块
**文件**：
- `R/hydro_pet.R`（重构自 PET.R）
- `R/hydro_balance.R`（重构自 WaterBalance.R）
- `R/hydro_analysis.R`（重构自 Hydro_obs.R）

**职责**：
- 蒸散发计算
- 水量平衡分析
- 水文观测处理

**关键函数**：
- `calc_pet_pm()` - Penman-Monteith PET
- `calc_melt_factor()` - 融雪因子
- `calc_water_balance()` - 水量平衡
- `plot_hydrograph()` - 水文过程线

#### 7. 可视化模块
**文件**：
- `R/plot_spatial.R`（重构自 plot.R, plotMap.R）
- `R/plot_timeseries.R`（新建）

**职责**：
- 空间数据可视化
- 时间序列可视化
- 对比图绘制

**关键函数**：
- `plot_mesh_2d()` - 2D 网格绘图（重构 map2d）
- `plot_timeseries()` - 时间序列绘图
- `compare_maps()` - 地图对比

#### 8. 主接口模块
**文件**：
- `R/interface_main.R`（重构自 autoBuildModel.R）

**职责**：
- 提供高级用户接口
- 协调各模块工作
- 自动格式转换

**关键函数**：
- `auto_build_model()` - 自动建模（重构 autoBuildModel）
- `quick_model()` - 快速建模（新增）

## 组件和接口

### 验证层接口设计

#### 参数验证接口

```r
# 验证正数
check_positive(x, name = "value", allow_zero = FALSE)
# 参数：
#   x: 数值
#   name: 参数名（用于错误消息）
#   allow_zero: 是否允许零
# 返回：invisible(TRUE) 或停止执行

# 验证文件存在
check_file_exists(path, name = "file", must_be_file = TRUE)
# 参数：
#   path: 文件路径
#   name: 参数名
#   must_be_file: 是否必须是文件（非目录）
# 返回：invisible(TRUE) 或停止执行

# 验证空间对象兼容性
check_spatial_compatible(x, y, check_crs = TRUE, check_overlap = FALSE,
                         name_x = "x", name_y = "y")
# 参数：
#   x, y: 空间对象（SpatRaster, SpatVector, sf）
#   check_crs: 是否检查 CRS
#   check_overlap: 是否检查范围重叠
# 返回：invisible(TRUE) 或停止执行

# CRS 兼容性检查
compatible_crs(crs1, crs2)
# 返回：logical
```

#### 直接使用 terra/sf 函数

rSHUD 不提供 terra/sf 的封装函数，直接使用原生 API：

```r
# 栅格操作 - 直接使用 terra
r <- terra::rast("dem.tif")
r_crop <- terra::crop(r, boundary)
slope <- terra::terrain(r_crop, "slope")

# 矢量操作 - 直接使用 sf
watershed <- sf::st_read("watershed.shp")
watershed_buf <- sf::st_buffer(watershed, dist = 1000)
area <- sf::st_area(watershed_buf)

# terra 和 sf 之间转换 - 使用原生函数
v <- terra::vect(sf_object)  # sf → SpatVector
sf_obj <- sf::st_as_sf(v)    # SpatVector → sf
```

### 核心函数接口设计

#### GIS 处理接口

```r
# 网格数据转栅格
mesh_to_raster(
  data,                    # 网格数据（向量或矩阵）
  mesh = NULL,             # 网格对象
  template = NULL,         # 栅格模板
  method = "idw",          # 插值方法：idw, linear, nearest
  resolution = NULL,       # 分辨率
  crs = NULL              # 坐标系统
)
# 返回：SpatRaster 对象

# 矢量转栅格
vector_to_raster(
  vector,                  # sf 或 SpatVector 对象
  template = NULL,         # 栅格模板
  field = 1,              # 栅格化字段
  resolution = NULL,       # 分辨率
  fun = "last"            # 聚合函数
)
# 返回：SpatRaster 对象
```

#### 网格生成接口

```r
# 生成三角网格
generate_mesh(
  boundary,                # 边界（sf polygon）
  rivers = NULL,          # 河流网络（sf lines）
  dem = NULL,             # 高程数据（SpatRaster）
  max_area = NULL,        # 最大三角形面积
  min_angle = 20,         # 最小角度
  extra_points = NULL     # 额外约束点
)
# 返回：list(
#   triangles = data.frame,  # 三角形定义
#   points = data.frame,     # 节点坐标
#   topology = list          # 拓扑关系
# )

# 计算网格属性
calc_mesh_attributes(
  mesh,                    # 网格对象
  soil = NULL,            # 土壤栅格
  geology = NULL,         # 地质栅格
  landcover = NULL,       # 土地覆盖栅格
  forcing = NULL,         # 强迫数据栅格
  ...                     # 其他属性栅格
)
# 返回：data.frame（属性表）
```

#### 河流处理接口

```r
# 构建河流网络
build_river_network(
  rivers,                  # 河流线（sf lines）
  dem,                    # 高程数据（SpatRaster）
  min_length = 100,       # 最小河段长度
  simplify = TRUE,        # 是否简化
  tolerance = 10          # 简化容差
)
# 返回：list(
#   network = sf,           # 河流网络
#   order = integer,        # 河流分级
#   downstream = integer,   # 下游关系
#   properties = data.frame # 河流属性
# )

# 计算河流属性（整合多个功能）
calc_river_properties(
  rivers,                  # 河流网络（sf lines）
  dem = NULL,             # 高程数据
  properties = c("order", "length", "slope", "downstream")
)
# 返回：sf 对象（带属性）
```

#### I/O 接口

```r
# 读取 SHUD 文件（统一接口）
read_shud(
  file,                    # 文件路径
  type = c("mesh", "river", "att", "para", "calib", "ic"),
  ...                     # 类型特定参数
)

# 写入 SHUD 文件（统一接口）
write_shud(
  data,                    # 数据对象
  file,                    # 文件路径
  type = c("mesh", "river", "att", "para", "calib", "ic"),
  ...                     # 类型特定参数
)

# 读取时间序列
read_timeseries(
  file,                    # 文件路径
  format = "tsd"          # 格式：tsd, csv
)
# 返回：xts 对象

# 读取模型输出
read_output(
  path,                    # 输出目录
  variables = NULL,       # 变量名
  time_range = NULL       # 时间范围
)
# 返回：list 或 xts
```

### 主接口设计

```r
# 自动建模主函数
auto_build_model(
  project_name,            # 项目名称
  domain,                  # 流域边界（sf polygon）
  rivers = NULL,          # 河流网络（sf lines）
  dem,                    # 高程数据（SpatRaster）
  soil = NULL,            # 土壤数据
  geology = NULL,         # 地质数据
  landcover = NULL,       # 土地覆盖数据
  forcing = NULL,         # 强迫数据
  output_dir = NULL,      # 输出目录
  mesh_options = list(),  # 网格生成选项
  river_options = list(), # 河流处理选项
  output_format = c("modern", "legacy"),  # 输出格式
  verbose = TRUE          # 是否显示进度
)
# 返回：list(
#   mesh = mesh_object,
#   river = river_object,
#   attributes = att_data,
#   files = file_paths
# )
```

## 数据模型

### 网格数据结构

```r
# SHUD.MESH S4 类（保持兼容，内部使用现代格式）
setClass("SHUD.MESH",
  slots = c(
    points = "data.frame",      # 节点坐标
    triangles = "data.frame",   # 三角形定义
    topology = "list",          # 拓扑关系
    attributes = "data.frame",  # 属性数据
    crs = "character",          # 坐标系统
    metadata = "list"           # 元数据
  )
)

# 内部使用 sf 格式存储
# - points: sf POINT
# - triangles: sf POLYGON
```

### 河流数据结构

```r
# SHUD.RIVER S4 类（保持兼容，内部使用现代格式）
setClass("SHUD.RIVER",
  slots = c(
    segments = "data.frame",    # 河段定义
    nodes = "data.frame",       # 节点
    topology = "list",          # 拓扑关系
    properties = "data.frame",  # 河流属性
    crs = "character",          # 坐标系统
    metadata = "list"           # 元数据
  )
)

# 内部使用 sf LINESTRING 格式
```

### 时间序列数据结构

```r
# 使用 xts 对象（保持不变）
# - 时间索引
# - 多列数据
# - 属性元数据
```

## 错误处理

### 错误处理策略

1. **参数验证**：在函数开始时验证所有参数
2. **类型检查**：使用 `inherits()` 检查对象类型
3. **值域检查**：验证数值参数在有效范围内
4. **兼容性检查**：验证空间对象的兼容性（CRS、范围等）
5. **清晰消息**：提供描述性错误消息

### 错误处理模板

```r
function_name <- function(x, y = NULL, tolerance = 0.01) {
  # 1. 必需参数检查
  if (missing(x)) {
    stop("参数 'x' 是必需的", call. = FALSE)
  }
  
  # 2. 类型检查
  if (!is_raster_like(x)) {
    stop(
      "参数 'x' 必须是栅格对象（SpatRaster 或 RasterLayer），",
      "但接收到 ", class(x)[1],
      call. = FALSE
    )
  }
  
  # 3. 值域检查
  if (tolerance <= 0) {
    stop(
      "参数 'tolerance' 必须为正数，但接收到 ", tolerance,
      call. = FALSE
    )
  }
  
  # 4. 兼容性检查
  if (!is.null(y)) {
    if (!compatible_crs(x, y)) {
      stop(
        "空间对象 CRS 不兼容：\n",
        "  x: ", get_crs(x), "\n",
        "  y: ", get_crs(y),
        call. = FALSE
      )
    }
  }
  
  # 5. 警告（非致命）
  if (tolerance > 1) {
    warning(
      "容差值较大 (", tolerance, ") 可能导致过度简化",
      call. = FALSE
    )
  }
  
  # 函数主体
  # ...
}
```

### 辅助验证函数

```r
# R/validators.R（新建）

# 检查正数
check_positive <- function(x, name = deparse(substitute(x))) {
  if (any(x <= 0)) {
    stop("'", name, "' 必须为正数", call. = FALSE)
  }
  invisible(TRUE)
}

# 检查文件存在
check_file_exists <- function(file) {
  if (!file.exists(file)) {
    stop("文件不存在: ", file, call. = FALSE)
  }
  invisible(TRUE)
}

# 检查空间兼容性
check_spatial_compatible <- function(x, y) {
  if (!compatible_geom(x, y)) {
    stop("空间对象不兼容", call. = FALSE)
  }
  invisible(TRUE)
}

# 检查 CRS 兼容性
compatible_crs <- function(x, y) {
  crs_x <- get_crs(x)
  crs_y <- get_crs(y)
  identical(crs_x, crs_y)
}
```

## 测试策略

### 测试层次

1. **单元测试**：测试单个函数
2. **集成测试**：测试模块间交互
3. **回归测试**：确保不破坏现有功能
4. **性能测试**：验证性能改进

### 测试组织

```
tests/
├── testthat/
│   ├── test-compat.R           # 兼容层测试
│   ├── test-gis-core.R         # GIS 核心功能
│   ├── test-mesh.R             # 网格生成
│   ├── test-river.R            # 河流处理
│   ├── test-io.R               # I/O 功能
│   ├── test-hydro.R            # 水文计算
│   ├── test-validation.R       # 参数验证
│   ├── test-integration.R      # 集成测试
│   └── helper-*.R              # 测试辅助函数
├── benchmark/
│   ├── bench-raster.R          # 栅格操作性能
│   ├── bench-vector.R          # 矢量操作性能
│   └── bench-overall.R         # 整体性能
└── fixtures/
    ├── test_data/              # 测试数据
    └── expected_results/       # 预期结果
```

### 测试示例

```r
# tests/testthat/test-compat.R

test_that("as_terra 正确转换 RasterLayer", {
  skip_if_not_installed("raster")
  skip_if_not_installed("terra")
  
  # 创建测试数据
  r_old <- raster::raster(ncol = 10, nrow = 10)
  raster::values(r_old) <- 1:100
  
  # 转换
  r_new <- as_terra(r_old)
  
  # 验证
  expect_s4_class(r_new, "SpatRaster")
  expect_equal(terra::ncell(r_new), 100)
  expect_equal(as.vector(terra::values(r_new)), 1:100)
})

test_that("mesh_to_raster 处理边界情况", {
  # 空输入
  expect_error(
    mesh_to_raster(NULL),
    "参数 'data' 是必需的"
  )
  
  # 无效数据类型
  expect_error(
    mesh_to_raster("invalid"),
    "必须是数值向量或矩阵"
  )
})

test_that("向后兼容性：旧函数仍然工作", {
  # 测试废弃函数
  expect_warning(
    result <- readmesh(test_file),
    "已废弃.*使用 read_mesh"
  )
  
  # 验证结果正确
  expect_s4_class(result, "SHUD.MESH")
})
```

### 性能测试

```r
# tests/benchmark/bench-raster.R

library(microbenchmark)

# 准备测试数据
n <- 1000
coords <- cbind(runif(n, 0, 100), runif(n, 0, 100))
values <- rnorm(n)

# 比较新旧实现
microbenchmark(
  legacy = {
    # 使用 raster 的旧实现
    r <- raster::rasterFromXYZ(cbind(coords, values))
  },
  modern = {
    # 使用 terra 的新实现
    r <- terra::rast(cbind(coords, values))
  },
  times = 100
)
```

## 实施计划

### 阶段 1：基础设施（第 1-2 周）

**目标**：建立兼容层和测试框架

**任务**：
1. 完善 `R/compat_terra_sf.R`
2. 创建 `R/validators.R`
3. 设置测试框架
4. 更新 DESCRIPTION 文件

**验收**：
- 兼容层函数全部实现并测试
- 测试框架可运行
- 包可以加载

### 阶段 2：核心 GIS 函数（第 3-4 周）

**目标**：迁移核心 GIS 处理函数

**任务**：
1. 重构 `mesh_to_raster()`
2. 重构 `vector_to_raster()`
3. 重构投影相关函数
4. 添加单元测试

**验收**：
- 所有 GIS 核心函数使用 terra/sf
- 测试覆盖率 > 80%
- 性能提升 > 20%

### 阶段 3：网格和河流模块（第 5-6 周）

**目标**：迁移网格生成和河流处理

**任务**：
1. 重构网格生成函数
2. 重构河流处理函数
3. 更新 S4 类定义
4. 添加集成测试

**验收**：
- 网格和河流模块完全迁移
- 与 GIS 模块集成正常
- demo 脚本可运行

### 阶段 4：I/O 和可视化（第 7 周）

**目标**：更新 I/O 和绘图函数

**任务**：
1. 重构 I/O 函数
2. 更新绘图函数
3. 确保输出格式兼容

**验收**：
- 可读写所有 SHUD 格式
- 绘图支持新格式
- 输出向后兼容

### 阶段 5：主接口和清理（第 8 周）

**目标**：更新主接口，清理冗余

**任务**：
1. 重构 `auto_build_model()`
2. 标记废弃函数
3. 清理冗余代码
4. 完善文档

**验收**：
- 主接口支持新旧格式
- 所有废弃函数标记
- 文档完整

### 阶段 6：测试和发布（第 9 周）

**目标**：全面测试和准备发布

**任务**：
1. 运行所有测试
2. 性能基准测试
3. 更新 NEWS.md
4. 准备发布材料

**验收**：
- R CMD check 通过
- 测试覆盖率 > 70%
- 文档完整
- 准备发布

## 性能优化

### 优化策略

1. **使用现代库**：terra 和 sf 本身就比旧库快
2. **向量化**：避免循环，使用向量化操作
3. **预分配**：为大对象预分配内存
4. **并行化**：对独立操作使用并行处理
5. **缓存**：缓存重复计算的结果

### 性能目标

| 操作类型 | 目标提升 | 测试数据规模 |
|---------|---------|------------|
| 栅格读写 | 2-5x | 1000x1000 |
| 矢量操作 | 2-10x | 10000 features |
| 网格生成 | 1.5-3x | 5000 triangles |
| 河流处理 | 2-5x | 1000 segments |
| 整体工作流 | 20-50% | 完整模型 |

## 文档策略

### 文档层次

1. **函数文档**：roxygen2 注释
2. **模块文档**：每个模块的概述
3. **用户指南**：vignettes
4. **迁移指南**：v2.x → v2.1
5. **开发文档**：贡献指南

### Vignettes 计划

1. **快速入门**：基本工作流
2. **模型构建**：详细的建模过程
3. **GIS 处理**：空间数据处理
4. **水文分析**：水文计算和分析
5. **迁移指南**：v2.1 → v2.2 升级指导

## 风险和缓解

### 主要风险

1. **破坏现有代码**
   - 缓解：兼容层 + 废弃警告 + 充分测试

2. **性能退化**
   - 缓解：性能基准测试 + 优化关键路径

3. **依赖问题**
   - 缓解：明确版本要求 + 可选依赖

4. **测试不足**
   - 缓解：多层次测试 + 高覆盖率目标

5. **文档滞后**
   - 缓解：同步更新 + 代码审查

## 成功标准

1. ✅ 所有核心函数直接使用 terra/sf（无封装层）
2. ✅ 完全移除 raster/sp/rgeos 依赖
3. ✅ 函数命名向后兼容（废弃别名保留）
4. ✅ 性能提升 ≥ 20%
5. ✅ 测试覆盖率 ≥ 70%
6. ✅ R CMD check 通过（0 错误 0 警告）
7. ✅ 所有 demo 使用 terra/sf 正常运行
8. ✅ 文档完整（所有导出函数）
9. ✅ 迁移指南完成（inst/MIGRATION_GUIDE.md）
10. ✅ 清晰的错误消息指导用户迁移
