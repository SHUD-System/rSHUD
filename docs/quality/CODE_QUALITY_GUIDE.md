# rSHUD 代码质量改进指南

## 1. 命名规范

### 当前问题
- 混合使用多种命名风格：`read.df()`, `readmesh()`, `sp.RiverOrder()`
- 不一致的前缀使用

### 改进方案

#### 函数命名
采用 **snake_case** 作为主要风格（符合 tidyverse 规范）：

```r
# 旧命名 → 新命名
read.df()        → read_df()
readmesh()       → read_mesh()
sp.RiverOrder()  → river_order() 或 sp_river_order()
MeshData2Raster()→ mesh_to_raster()
```

#### 前缀规范
- `read_*` - 读取函数
- `write_*` - 写入函数
- `calc_*` - 计算函数
- `get_*` - 获取/提取函数
- `set_*` - 设置函数
- `check_*` - 验证函数
- `plot_*` - 绘图函数

#### 变量命名
```r
# 好的命名
n_cells <- 100
river_length <- 1500
dem_resolution <- 30

# 避免的命名
nx <- 100  # 不清晰
rl <- 1500  # 缩写不明确
```

### 迁移策略
```r
# 保留旧函数名，添加废弃警告
read.df <- function(...) {
  .Deprecated("read_df")
  read_df(...)
}

# 新函数
read_df <- function(file, ...) {
  # 实现
}
```

## 2. 参数验证和错误处理

### 当前问题
- 缺少参数类型检查
- 错误信息不够清晰
- 没有输入验证

### 改进方案

#### 参数验证模板
```r
function_name <- function(x, y = NULL, tolerance = 0.01) {
  # 1. 必需参数检查
  if (missing(x)) {
    stop("Argument 'x' is required", call. = FALSE)
  }
  
  # 2. 类型检查
  if (!inherits(x, c("SpatRaster", "RasterLayer"))) {
    stop("'x' must be a raster object (SpatRaster or RasterLayer)", 
         call. = FALSE)
  }
  
  # 3. 值域检查
  if (tolerance <= 0) {
    stop("'tolerance' must be positive", call. = FALSE)
  }
  
  # 4. 逻辑检查
  if (!is.null(y) && !compatible(x, y)) {
    stop("'x' and 'y' must have compatible dimensions", call. = FALSE)
  }
  
  # 5. 警告（非致命问题）
  if (tolerance > 1) {
    warning("Large tolerance value may result in over-simplification",
            call. = FALSE)
  }
  
  # 函数主体
  # ...
}
```

#### 错误信息规范
```r
# ❌ 不好的错误信息
stop("Error!")
stop("Invalid input")

# ✅ 好的错误信息
stop("Invalid raster resolution: expected positive number, got ", res)
stop("File not found: ", file_path)
stop("Incompatible CRS: x has '", crs(x), "', y has '", crs(y), "'")
```

#### 辅助验证函数
```r
# 创建 R/validators.R
check_positive <- function(x, name = deparse(substitute(x))) {
  if (any(x <= 0)) {
    stop("'", name, "' must be positive", call. = FALSE)
  }
  invisible(TRUE)
}

check_file_exists <- function(file) {
  if (!file.exists(file)) {
    stop("File not found: ", file, call. = FALSE)
  }
  invisible(TRUE)
}

check_spatial_compatible <- function(x, y) {
  if (!terra::compareGeom(x, y, stopOnError = FALSE)) {
    stop("Spatial objects are not compatible", call. = FALSE)
  }
  invisible(TRUE)
}
```

## 3. 函数文档规范

### 当前问题
- 部分函数缺少文档
- 文档不完整（缺少参数说明、返回值）
- 缺少示例

### 改进方案

#### 完整文档模板
```r
#' Calculate river network properties
#'
#' This function calculates various properties of a river network including
#' stream order, length, slope, and downstream relationships.
#'
#' @param river_lines An sf or SpatialLines object representing the river network
#' @param dem A SpatRaster or RasterLayer object with elevation data
#' @param min_length Minimum river segment length in meters. Segments shorter
#'   than this will be merged with adjacent segments. Default is 100.
#' @param simplify Logical. Should the river network be simplified before
#'   processing? Default is TRUE.
#' @param tolerance Simplification tolerance in map units. Only used if
#'   \code{simplify = TRUE}. Default is 10.
#'
#' @return A list with the following components:
#'   \item{river}{An sf object with river segments and attributes}
#'   \item{order}{Integer vector of stream orders (Strahler)}
#'   \item{length}{Numeric vector of segment lengths in meters}
#'   \item{slope}{Numeric vector of segment slopes (m/m)}
#'   \item{downstream}{Integer vector of downstream segment IDs}
#'
#' @details
#' The function uses the Strahler method for stream ordering. River slopes
#' are calculated from the DEM using bilinear interpolation at segment
#' endpoints. Segments with negative slopes are automatically corrected.
#'
#' @section Migration Note:
#' As of version 3.0, this function uses \code{sf} and \code{terra} internally.
#' Old \code{sp} and \code{raster} objects are automatically converted.
#'
#' @seealso
#' \code{\link{river_order}} for stream ordering only,
#' \code{\link{river_slope}} for slope calculation only
#'
#' @export
#' @examples
#' \dontrun{
#' # Load example data
#' library(terra)
#' library(sf)
#' 
#' dem <- rast(system.file("extdata/dem.tif", package = "rSHUD"))
#' rivers <- st_read(system.file("extdata/rivers.shp", package = "rSHUD"))
#' 
#' # Calculate river properties
#' result <- calc_river_properties(rivers, dem)
#' 
#' # Plot stream order
#' plot(result$river["order"], main = "Stream Order")
#' 
#' # Summary statistics
#' summary(result$slope)
#' }
#'
#' @references
#' Strahler, A. N. (1957). Quantitative analysis of watershed geomorphology.
#' Transactions, American Geophysical Union, 38(6), 913-920.
calc_river_properties <- function(river_lines, dem, 
                                   min_length = 100,
                                   simplify = TRUE,
                                   tolerance = 10) {
  # 实现
}
```

#### 文档检查清单
- [ ] 标题（一句话描述）
- [ ] 详细描述（多段落）
- [ ] 所有参数说明（@param）
- [ ] 返回值说明（@return）
- [ ] 详细说明（@details）
- [ ] 示例代码（@examples）
- [ ] 相关函数（@seealso）
- [ ] 导出标记（@export）
- [ ] 参考文献（@references，如适用）

## 4. 代码组织和模块化

### 当前问题
- 单个文件过长（如 `Func_GIS.R` 超过 700 行）
- 函数职责不清晰
- 重复代码

### 改进方案

#### 文件拆分原则
每个文件应该：
- 包含相关功能的函数（单一职责）
- 不超过 500 行
- 有清晰的文件头注释

```r
# R/river_network.R
# River network processing functions
# 
# This file contains functions for processing river networks including:
# - Stream ordering
# - Slope calculation
# - Network topology
# - Segment merging

# R/river_order.R - 河流分级
# R/river_slope.R - 坡度计算
# R/river_topology.R - 拓扑关系
```

#### 函数拆分原则
```r
# ❌ 不好：一个函数做太多事
process_river <- function(river, dem, ...) {
  # 100+ 行代码
  # 做了：简化、分级、计算坡度、拓扑分析...
}

# ✅ 好：拆分成多个小函数
process_river <- function(river, dem, ...) {
  river <- simplify_river(river, ...)
  order <- calc_stream_order(river)
  slope <- calc_river_slope(river, dem)
  topo <- build_topology(river)
  
  list(river = river, order = order, slope = slope, topology = topo)
}

simplify_river <- function(river, tolerance) {
  # 20 行代码
}

calc_stream_order <- function(river) {
  # 30 行代码
}

calc_river_slope <- function(river, dem) {
  # 40 行代码
}

build_topology <- function(river) {
  # 50 行代码
}
```

#### 避免重复代码
```r
# ❌ 不好：重复的坐标提取代码
func1 <- function(x) {
  if (inherits(x, "sf")) {
    coords <- sf::st_coordinates(x)
  } else if (inherits(x, "SpatVector")) {
    coords <- terra::crds(x)
  }
  # ...
}

func2 <- function(x) {
  if (inherits(x, "sf")) {
    coords <- sf::st_coordinates(x)
  } else if (inherits(x, "SpatVector")) {
    coords <- terra::crds(x)
  }
  # ...
}

# ✅ 好：提取为辅助函数
extract_coords <- function(x) {
  if (inherits(x, "sf")) {
    return(sf::st_coordinates(x))
  } else if (inherits(x, "SpatVector")) {
    return(terra::crds(x))
  }
  stop("Unsupported object type")
}

func1 <- function(x) {
  coords <- extract_coords(x)
  # ...
}

func2 <- function(x) {
  coords <- extract_coords(x)
  # ...
}
```

## 5. 性能优化

### 向量化操作
```r
# ❌ 慢：循环
result <- numeric(n)
for (i in 1:n) {
  result[i] <- sqrt(x[i]^2 + y[i]^2)
}

# ✅ 快：向量化
result <- sqrt(x^2 + y^2)
```

### 预分配内存
```r
# ❌ 慢：动态增长
result <- c()
for (i in 1:n) {
  result <- c(result, compute(i))
}

# ✅ 快：预分配
result <- vector("numeric", n)
for (i in 1:n) {
  result[i] <- compute(i)
}
```

### 使用合适的数据结构
```r
# ❌ 慢：data.frame 用于大量行操作
df <- data.frame(x = 1:1e6, y = 1:1e6)
df$z <- df$x + df$y

# ✅ 快：matrix 用于数值计算
mat <- matrix(c(1:1e6, 1:1e6), ncol = 2)
z <- mat[, 1] + mat[, 2]
```

## 6. 测试规范

### 单元测试模板
```r
# tests/testthat/test-river-order.R

test_that("river_order calculates correct Strahler order", {
  # 准备测试数据
  river <- create_test_river()
  
  # 执行函数
  result <- river_order(river)
  
  # 断言
  expect_type(result, "integer")
  expect_length(result, nrow(river))
  expect_true(all(result >= 1))
  expect_equal(max(result), 3)  # 已知测试数据的最大级别
})

test_that("river_order handles edge cases", {
  # 空输入
  expect_error(river_order(NULL), "required")
  
  # 单条河流
  single_river <- create_single_river()
  expect_equal(river_order(single_river), 1)
  
  # 无连接的河流
  disconnected <- create_disconnected_rivers()
  result <- river_order(disconnected)
  expect_true(all(result == 1))
})

test_that("river_order works with different input formats", {
  river_sf <- create_test_river_sf()
  river_sp <- create_test_river_sp()
  
  result_sf <- river_order(river_sf)
  result_sp <- river_order(river_sp)
  
  expect_equal(result_sf, result_sp)
})
```

## 7. 代码审查清单

在提交代码前检查：

### 功能性
- [ ] 函数按预期工作
- [ ] 处理边界情况
- [ ] 错误处理适当
- [ ] 参数验证完整

### 代码质量
- [ ] 命名清晰一致
- [ ] 没有重复代码
- [ ] 函数职责单一
- [ ] 代码可读性好

### 文档
- [ ] 函数文档完整
- [ ] 参数说明清楚
- [ ] 有工作示例
- [ ] 更新了 NEWS.md

### 测试
- [ ] 有单元测试
- [ ] 测试覆盖主要路径
- [ ] 测试边界情况
- [ ] 所有测试通过

### 性能
- [ ] 没有明显的性能问题
- [ ] 使用了向量化操作
- [ ] 内存使用合理

### 兼容性
- [ ] 向后兼容（或有迁移路径）
- [ ] 支持多种输入格式
- [ ] 依赖包版本合理

## 8. 持续改进

### 代码度量
定期检查：
- 函数长度（目标：< 50 行）
- 文件长度（目标：< 500 行）
- 圈复杂度（目标：< 10）
- 测试覆盖率（目标：> 80%）

### 工具
```r
# 代码风格检查
lintr::lint_package()

# 测试覆盖率
covr::package_coverage()

# 代码复杂度
cyclocomp::cyclocomp_package()
```

### 重构优先级
1. 🔴 高：核心函数、频繁使用的函数
2. 🟡 中：辅助函数、中等使用频率
3. 🟢 低：很少使用的函数、即将废弃的函数
