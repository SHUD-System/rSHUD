# rSHUD --- SHUD建模系统工具箱

[![R-CMD-check](https://github.com/SHUD-System/rSHUD/workflows/R-CMD-check/badge.svg)](https://github.com/SHUD-System/rSHUD/actions)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R](https://img.shields.io/badge/R-≥3.5.0-blue.svg)](https://www.r-project.org/)

**作者**: Lele Shu (shulele@lzb.ac.cn)  
**官网**: [www.shud.xyz](https://www.shud.xyz/)  
**GitHub**: [SHUD-System/rSHUD](https://github.com/SHUD-System/rSHUD)

---

## 📖 项目简介

**rSHUD** 是为SHUD（Simulator for Hydrological Unstructured Domain，水文非结构化域模拟器）水文模型系统设计的R语言工具箱。

### SHUD模型特点
SHUD是一个多过程、多尺度的水文模型，使用半离散有限体积方法将主要水文过程完全耦合。该包可与AutoSHUD项目配合使用，实现建模域的自动构建。

## 📚 项目文档索引

流程与治理类文档已整理到 `docs/` 目录：
- [文档索引](docs/README.md)

---

## 🎯 主要功能

### 1. 数据转换与预处理
- 将地理空间数据转换为SHUD格式
- 处理栅格和矢量数据
- 构建非结构化三角网格域

### 2. 模型文件管理
- 读写SHUD输入文件
- 读取SHUD输出文件
- 生成校准参数集

### 3. 水文分析
- 水文数据时间序列分析
- 二维和三维可视化绘图
- 水文图分析

### 4. GIS功能
- 将非结构化数据转换为空间数据（Shapefile或Raster）
- 流域划分与河流网络处理
- 空间分析和可视化

---

## 🚀 安装说明

### 系统依赖（Ubuntu平台）
```bash
sudo apt -y install gdal-bin libgdal-dev gcc g++ gfortran
sudo apt -y install r-cran-systemfonts r-cran-textshaping
```

### 安装依赖包
```r
# 建议先安装依赖包
libs = c('reshape2','xts','hydroGOF','zoo','RTriangle','proj4','gstat',
         'abind','lubridate','interp','geometry','testthat', 'rmarkdown', 
         'ncdf4', 'blogdown', 'doParallel', 'knitr', 'rmarkdown', 'deldir',
         'ncdf4', 'devtools')

for(i in 1:length(libs)){
  if(!require(libs[i], character.only = TRUE)){
    install.packages(libs[i], dependencies = TRUE, INSTALL_opts = '--no-lock')
  }
}

# 安装RTriangle包
if(!require(devtools)){
  install.packages("devtools", dependencies = TRUE, INSTALL_opts = '--no-lock')
}
devtools::install_github("shulele/RTriangle", subdir="pkg")
```

### 安装rSHUD
```r
if(!require(devtools)){ install.packages("devtools") }
devtools::install_github("SHUD-System/rSHUD")
```

---

## 📚 使用示例

### 基本使用
```r
library(rSHUD)

# 加载示例数据
data(sac)
indata = sac

# 自动构建模型
sp.forc = indata[['forc']]
forc.fns = paste0(sp.forc@data[, 'NLDAS_ID'], '.csv')
pm = autoBuildModel(sac, forcfiles = forc.fns, outdir='./output')

# 转换为sf对象
spm = mesh_to_sf(pm)
```

### 查看可用示例
```r
# 查看所有示例脚本
list.files(system.file("demo", package = "rSHUD"))

# 运行特定示例
demo(demo_autoBuild, package = "rSHUD")
```

---

## 📁 项目结构

```
rSHUD/
├── R/                    # R源代码
├── man/                  # 函数帮助文档
├── demo/                 # 示例脚本
├── data/                 # 示例数据集
├── src/                  # C++源代码
├── experiments/          # 实验代码
└── output/               # 输出目录
```

---

## 🔧 主要函数分类

### 模型构建
- `autoBuildModel()` - 自动构建SHUD模型
- `readmesh()`, `readriv()` - 读取网格和河流数据
- `write_mesh()`, `write_river()` - 写入网格和河流数据

### GIS功能
- `watershedDelineation()` - 流域划分
- `sp2raster()` - 空间数据转栅格
- `writeshape()` - 输出Shapefile

### 水文计算
- `PET_PM()` - 彭曼-蒙蒂斯蒸散发计算
- `MeltFactor()` - 融雪因子
- `hydrograph()` - 水文图分析

---

## 📖 文档资源

- **函数帮助**: 使用 `?function_name` 查看详细帮助
- **示例代码**: 查看 `demo/` 目录中的示例脚本
- **在线文档**: [www.shud.xyz](https://www.shud.xyz/)
- **英文版本**: [README.md](README.md)

---

## 🤝 贡献指南

欢迎提交Issue和Pull Request！请确保：
1. 代码符合R包开发规范
2. 新函数包含完整的文档和测试
3. 遵循项目的编码风格

详细贡献指南请参见 [CONTRIBUTING.md](CONTRIBUTING.md)。

---

## 📄 许可证

本项目采用MIT许可证 - 详见 [LICENSE](LICENSE) 文件

---

## 📞 联系方式

- **邮箱**: shulele@lzb.ac.cn
- **项目主页**: [www.shud.xyz](https://www.shud.xyz/)
- **GitHub**: [SHUD-System/rSHUD](https://github.com/SHUD-System/rSHUD)

---

*如果这个项目对您有帮助，请给我们一个⭐️星标！*
