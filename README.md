# 高效 TDOA 算法实现

基于 C 语言的高效 TDOA (Time Difference of Arrival) 定位算法实现。

## 特性

- **Fang 算法**: 闭式解算法，适用于 2D 定位
- **Chan 算法**: 梯度下降法，适用于 2D/3D 定位
- **Taylor 级数法**: 高斯-牛顿迭代法，需要好的初始猜测
- **C99/C11 标准**: 符合现代 C 语言规范
- **正确内存管理**: malloc/free，无内存泄漏

## 构建

```bash
make          # 编译
make test     # 运行演示
make clean    # 清理
```

## 使用示例

```c
#include "tdoa.h"

tdoa_solver_t solver;
tdoa_solver_init(&solver);

// 添加接收器
tdoa_point3d_t r1 = {0.0, 0.0, 0.0};
tdoa_point3d_t r2 = {100.0, 0.0, 0.0};
tdoa_point3d_t r3 = {50.0, 86.6, 0.0};

tdoa_add_receiver(&solver, &r1, 0.0);
tdoa_add_receiver(&solver, &r2, 0.0);
tdoa_add_receiver(&solver, &r3, 0.0);

// 生成测量数据
tdoa_point3d_t true_pos = {50.0, 30.0, 0.0};
tdoa_generate_measurements(&solver, &true_pos, 1e-9);

// 求解位置
tdoa_point3d_t result;
tdoa_solve_chan(&solver, &result);

tdoa_solver_destroy(&solver);
```

## API

### 初始化

- `tdoa_solver_init()` - 初始化求解器
- `tdoa_solver_configure()` - 配置参数
- `tdoa_solver_destroy()` - 释放资源

### 添加数据

- `tdoa_add_receiver()` - 添加接收器
- `tdoa_add_measurement()` - 添加 TDOA 测量

### 求解算法

- `tdoa_solve_fang()` - Fang 闭式解
- `tdoa_solve_chan()` - Chan 梯度下降
- `tdoa_solve_taylor()` - Taylor 级数迭代

### 工具

- `tdoa_generate_measurements()` - 生成模拟测量数据
- `tdoa_get_speed_of_light()` - 光速
- `tdoa_get_speed_of_sound()` - 声速

## 测试结果

| 测试 | 算法 | 误差 |
|------|------|------|
| 2D (3接收器) | Chan | ~0.33m |
| 3D (4接收器) | Chan/Taylor | ~0.22m |
| Monte Carlo | Taylor | ~0.3m (1ns噪声) |

## 许可

MIT License
