# 贡献指南

感谢您对rSHUD项目的关注！我们欢迎所有形式的贡献，包括但不限于代码贡献、文档改进、问题报告和功能建议。

---

## 🚀 如何贡献

### 1. 报告问题
如果您发现了bug或有功能建议，请通过GitHub Issues提交：

- 使用清晰、描述性的标题
- 提供详细的复现步骤
- 包含系统信息和R版本
- 如果可能，提供可重现的代码示例

### 2. 提交代码
我们欢迎代码贡献！请遵循以下步骤：

#### 2.1 Fork项目
1. 在GitHub上fork rSHUD项目
2. 克隆您的fork到本地
3. 创建新的功能分支

```bash
git clone https://github.com/YOUR_USERNAME/rSHUD.git
cd rSHUD
git checkout -b feature/your-feature-name
```

#### 2.2 开发规范
- **代码风格**: 遵循R语言编码规范
- **文档**: 所有新函数必须包含完整的Roxygen2文档
- **测试**: 为新功能添加相应的测试用例
- **中英文支持**: 所有用户可见的文本都应提供中英文版本

#### 2.3 提交代码
```bash
git add .
git commit -m "feat: add new function for watershed analysis"
git push origin feature/your-feature-name
```

### 3. 文档贡献
我们特别欢迎文档改进：

- 完善函数帮助文档
- 添加使用示例和教程
- 改进README和项目说明
- 翻译和本地化

---

## 📋 代码规范

### R代码规范
- 使用4个空格缩进
- 函数名使用驼峰命名法（如：`autoBuildModel`）
- 变量名使用点分隔（如：`data.input`）
- 每行代码不超过80个字符

### 文档规范
```r
#' 函数标题
#' 
#' 函数描述
#' 
#' @param param1 参数1描述
#' @param param2 参数2描述
#' @return 返回值描述
#' @export
#' @examples
#' # 示例代码
#' function_name(param1, param2)
```

---

## 🧪 测试要求

### 单元测试
- 所有新功能必须包含测试用例
- 使用`testthat`包编写测试
- 测试覆盖率应达到80%以上

### 测试运行
```r
# 运行所有测试
devtools::test()

# 检查代码覆盖率
devtools::test_coverage()
```

---

## 📝 提交信息规范

我们使用约定式提交信息格式：

- `feat:` 新功能
- `fix:` 修复bug
- `docs:` 文档更新
- `style:` 代码格式调整
- `refactor:` 代码重构
- `test:` 测试相关
- `chore:` 构建过程或辅助工具的变动

示例：
```
feat: add watershed delineation function with bilingual support
fix: resolve memory leak in large dataset processing
docs: update README with installation instructions
```

---

## 🔄 工作流程

### 1. 讨论功能
在开始编码之前，建议先在Issues中讨论新功能：

- 描述功能需求
- 讨论实现方案
- 获得社区反馈

### 2. 开发分支
- 从`main`分支创建功能分支
- 定期同步主分支的更新
- 保持分支的整洁和专注

### 3. 代码审查
- 所有代码必须通过审查
- 审查者会检查代码质量、测试覆盖率和文档完整性
- 根据反馈进行必要的修改

---

## 🎯 优先贡献领域

我们特别欢迎以下领域的贡献：

### 高优先级
- 性能优化和内存管理
- 错误处理和用户友好的错误信息
- 单元测试和集成测试
- 中英文文档完善

### 中优先级
- 新算法和方法的实现
- 可视化功能的增强
- 示例代码和教程
- 包依赖的现代化

---

## 📞 获取帮助

如果您在贡献过程中遇到问题：

- 查看现有的Issues和Pull Requests
- 在Issues中提问
- 联系维护者：shulele@lzb.ac.cn

---

## 🙏 致谢

感谢所有为rSHUD项目做出贡献的开发者！您的贡献使这个项目变得更好。

---

*让我们一起构建更好的水文建模工具！*

---

**English Version**: [CONTRIBUTING.md](CONTRIBUTING.md)
