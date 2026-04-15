# rSHUD 升级实施指南

## 开始之前

### 1. 创建开发分支
```bash
git checkout -b dev-v3.0
```

### 2. 更新 DESCRIPTION 文件

将以下内容更新到 `DESCRIPTION`：

```r
Depends: R (>= 4.0.0)
Imports: 
    Rcpp,
    terra (>= 1.7-0),
    sf (>= 1.0-0),
    reshape2,
    ggplot2,
    gridExtra,
    grid,
    fields,
    xts, 
    hydroGOF,
    zoo,
    RTriangle,
    gstat,
    abind,
    utils,
    lubridate,
    doParallel,
    geometry,
    methods
Suggests: 
    testthat (>= 3.0.0),
    knitr,
    rmarkdown,
    deldir,
    interp,
    whitebox,
    ncdf4,
    raster,
    sp,
    rgeos
```

**关键变更**：
- 将 `terra` 和 `sf` 从 Suggests 移到 Imports
- 将 `raster`, `sp`, `rgeos` 移到 Suggests（向后兼容）
- 提高 R 版本要求到 4.0.0

### 3. 安装新依赖

```r
install.packages(c("terra", "sf"))
```

## 迁移步骤

### 步骤 1：创建兼容层 ✅

已完成：`R/compat_terra_sf.R`

这个文件提供了统一的接口函数：
- `as_terra()` - 转换为 terra 格式
- `as_sf()` - 转换为 sf 格式
- `get_extent()` - 获取空间范围
- `get_crs()` / `set_crs()` - CRS 操作
- `simplify_geom()` - 几何简化
- `calc_area()` - 面积计算
- `extract_coords()` - 坐标提取

### 步骤 2：迁移核心 GIS 函数

#### 2.1 更新 `MeshData2Raster()` 函数

**位置**：`R/Func_GIS.R`

**当前实现**（使用 raster）：
```r
MeshData2Raster <- function(x, rmask = NULL, ...) {
  # ... 省略 ...
  r = raster::raster(dat)
  # ... 省略 ...
}
```

**新实现**（使用 terra）：
```r
MeshData2Raster <- function(x, rmask = NULL, method = 'idw', stack = FALSE) {
  # 参数验证
  if (!is.null(rmask) && !is_raster_like(rmask)) {
    stop("rmask must be a raster object (SpatRaster or RasterLayer)")
  }
  
  # 转换为 terra 格式
  if (!is.null(rmask)) {
    rmask <- as_terra(rmask)
  }
  
  # ... 核心逻辑 ...
  
  # 使用 terra 函数
  r <- terra::rast(dat)
  
  return(r)
}
```

#### 2.2 更新 `sp2raster()` 函数

**新实现**：
```r
sp2raster <- function(sp, mask = NULL, resolution = NULL, field = 1) {
  # 转换为 sf 格式
  if (!inherits(sp, "sf")) {
    sp <- as_sf(sp)
  }
  
  # 创建栅格模板
  if (is.null(mask)) {
    ext <- sf::st_bbox(sp)
    if (is.null(resolution)) {
      resolution <- min(ext["xmax"] - ext["xmin"], 
                       ext["ymax"] - ext["ymin"]) / 100
    }
    r <- terra::rast(ext, resolution = resolution)
  } else {
    r <- as_terra(mask)
  }
  
  # 栅格化
  r <- terra::rasterize(terra::vect(sp), r, field = field)
  
  return(r)
}
```

### 步骤 3：迁移河流处理函数

#### 3.1 更新 `sp.RiverOrder()` 函数

**位置**：`R/GIS_RiverProcess.R`

**新实现**：
```r
sp.RiverOrder <- function(sp, coord = NULL) {
  msg <- 'sp.RiverOrder::'
  
  # 转换为 sf 格式
  if (!inherits(sp, "sf")) {
    sp <- as_sf(sp)
  }
  
  # 简化几何
  ext <- sf::st_bbox(sp)
  tol <- (ext["xmax"] - ext["xmin"]) * 0.01
  sp <- sf::st_simplify(sp, dTolerance = tol)
  
  # 提取坐标
  if (is.null(coord)) {
    coord <- extract_coords(sp, unique = TRUE)
  }
  
  # ... 其余逻辑保持不变 ...
}
```

### 步骤 4：更新主接口函数

#### 4.1 更新 `autoBuildModel()` 函数

**位置**：`R/autoBuildModel.R`

**策略**：
1. 自动检测输入格式
2. 内部统一使用 terra/sf
3. 输出可选择格式

```r
autoBuildModel <- function(
  prjname,
  domain,
  riv = NULL,
  dem,
  # ... 其他参数 ...
  output_format = c("terra", "raster")  # 新参数
) {
  output_format <- match.arg(output_format)
  
  # 自动转换输入
  if (!inherits(domain, "sf")) {
    message("Converting domain to sf format...")
    domain <- as_sf(domain)
  }
  
  if (!is.null(riv) && !inherits(riv, "sf")) {
    message("Converting river to sf format...")
    riv <- as_sf(riv)
  }
  
  if (!inherits(dem, "SpatRaster")) {
    message("Converting DEM to terra format...")
    dem <- as_terra(dem)
  }
  
  # ... 核心逻辑 ...
  
  # 根据需要转换输出格式
  if (output_format == "raster" && requireNamespace("raster", quietly = TRUE)) {
    # 转换回旧格式（向后兼容）
    result$dem <- raster::raster(result$dem)
    # ... 其他转换 ...
  }
  
  return(result)
}
```

### 步骤 5：更新文档

#### 5.1 添加迁移指南到 README

在 `README.md` 中添加：

```markdown
## Version 3.0 Migration Guide

rSHUD v3.0 has migrated from `raster`/`sp`/`rgeos` to `terra`/`sf` for better performance and modern spatial data handling.

### For New Users

Simply install and use:
```r
install.packages("rSHUD")
library(rSHUD)
```

### For Existing Users

Your old code will mostly work, but we recommend updating:

**Old code (v2.x)**:
```r
library(raster)
library(sp)
dem <- raster("dem.tif")
wbd <- readOGR("watershed.shp")
```

**New code (v3.x)**:
```r
library(terra)
library(sf)
dem <- rast("dem.tif")
wbd <- st_read("watershed.shp")
```

The main functions (`autoBuildModel`, `shud.triangle`, etc.) automatically detect and convert input formats.

### Compatibility

- Old `raster`/`sp` objects are automatically converted
- Use `output_format = "raster"` in main functions to get old format outputs
- Deprecated functions will show warnings but still work
```

#### 5.2 更新函数文档

为每个迁移的函数添加：

```r
#' @section Migration Note:
#' As of version 3.0, this function uses \code{terra} and \code{sf} internally.
#' Old \code{raster} and \code{sp} objects are automatically converted.
#' To get output in old format, use \code{output_format = "raster"}.
```

### 步骤 6：添加测试

创建 `tests/testthat/test-compat.R`：

```r
test_that("Compatibility layer works", {
  skip_if_not_installed("terra")
  skip_if_not_installed("sf")
  skip_if_not_installed("raster")
  skip_if_not_installed("sp")
  
  # Test raster conversion
  r_old <- raster::raster(ncol = 10, nrow = 10)
  r_new <- as_terra(r_old)
  expect_s4_class(r_new, "SpatRaster")
  
  # Test sp conversion
  coords <- cbind(1:10, 1:10)
  sp_old <- sp::SpatialPoints(coords)
  sf_new <- as_sf(sp_old)
  expect_s3_class(sf_new, "sf")
  
  # Test extent extraction
  ext <- get_extent(r_old)
  expect_length(ext, 4)
  expect_type(ext, "double")
})
```

### 步骤 7：性能测试

创建 `tests/benchmark.R`：

```r
# 比较新旧实现的性能
library(microbenchmark)

# 测试数据
n <- 1000
coords <- cbind(runif(n, 0, 100), runif(n, 0, 100))
values <- rnorm(n)

# 旧实现 vs 新实现
microbenchmark(
  old = {
    # 使用 raster 的旧代码
  },
  new = {
    # 使用 terra 的新代码
  },
  times = 10
)
```

## 测试清单

在每个迁移步骤后，运行以下测试：

```r
# 1. 加载包
devtools::load_all()

# 2. 运行测试
devtools::test()

# 3. 检查文档
devtools::document()
devtools::check()

# 4. 运行示例
demo(demo_autoBuild, package = "rSHUD")
```

## 常见问题

### Q1: 如何处理依赖旧包的用户代码？

A: 提供兼容层和自动转换。用户可以继续传入旧格式对象，函数内部自动转换。

### Q2: 性能会有提升吗？

A: 是的，terra 通常比 raster 快 2-10 倍，特别是大数据集。

### Q3: 需要重写所有代码吗？

A: 不需要。核心函数内部迁移，用户接口保持兼容。

### Q4: 如何处理 Rcpp 代码？

A: Rcpp 代码不需要改动，继续使用。只需确保 R 层面的数据转换正确。

## 发布计划

1. **Alpha 版本** (内部测试)
   - 完成核心函数迁移
   - 基本测试通过

2. **Beta 版本** (公开测试)
   - 完成所有函数迁移
   - 文档更新完成
   - 征求用户反馈

3. **RC 版本** (候选发布)
   - 修复 Beta 版本问题
   - 性能优化
   - 最终测试

4. **v3.0.0 正式版**
   - 发布到 CRAN
   - 更新网站文档
   - 发布公告

## 回滚计划

如果遇到严重问题：

1. 保留 v2.x 分支
2. 标记 v3.0 为 experimental
3. 继续维护 v2.x 直到 v3.0 稳定

## 资源

- [terra 文档](https://rspatial.github.io/terra/)
- [sf 文档](https://r-spatial.github.io/sf/)
- [从 raster 迁移到 terra](https://rspatial.org/spatial/8-rastermanip.html)
- [从 sp 迁移到 sf](https://r-spatial.github.io/sf/articles/sf1.html)
