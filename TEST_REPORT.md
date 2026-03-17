# 🔍 TDOA Demo 测试报告

**项目**: tdoa-demo (bwandmff/tdoa-demo)  
**日期**: 2026-03-17  
**测试工程师**: 严苛测试模式 🔧

---

## 📊 测试执行结果

| 测试 | 状态 | 备注 |
|------|------|------|
| 编译 | ✅ 通过 | gcc + Wall -Wextra -Wpedantic |
| 运行 | ✅ 通过 | 所有测试用例执行完成 |

---

## 🐛 发现的问题

### 🔴 严重问题 (Critical)

#### 1. Monte Carlo 测试误差异常 (100米误差)

- **位置**: `main.c` Test 3
- **现象**: 1ns噪声下平均误差100米（理论应~0.3米）
- **根因**: `tdoa_solver_reset()` 释放了 `measurements` 内存并将 `num_receivers` 设为0，但循环中只调用了 `tdoa_add_receiver` 没有重新初始化求解器
- **代码位置**: `main.c:164-178`

```c
// 问题：reset后num_receivers=0，但add_receiver会累加
tdoa_solver_reset(&solver);
tdoa_add_receiver(&solver, &r1, 0.0);  // num_receivers从1开始，而不是0！
```

---

#### 2. 声学TDOA使用错误的速度常量

- **位置**: `main.c` Test 4 + `tdoa.c:tdoa_generate_measurements()`
- **现象**: 声速测试误差1666米
- **根因**: 测量生成时使用 `SPEED_OF_LIGHT` (299,792,458 m/s)，而不是声速 (343 m/s)
- **代码位置**: `tdoa.c:554` `const double c = SPEED_OF_LIGHT;`

```c
// 应该根据应用场景选择速度！
const double c = SPEED_OF_LIGHT;  // 应用于声学TDOA是错误的
```

---

#### 3. 内存释放后使用 (静态分析警告)

- **位置**: `tdoa.c:193`
- **警告**: `use of memory after it is freed`
- **问题**: `safe_realloc` 内部先 `free(ptr)` 再分配，但返回前使用了传入的 `ptr`

```c
static void *safe_realloc(void *ptr, size_t size) {
    void *new_ptr = realloc(ptr, size);
    if (new_ptr == NULL) {
        free(ptr);  // 如果realloc失败，这里free了ptr
        // 但调用者可能还在使用ptr！
    }
    return new_ptr;
}
```

---

### 🟠 中等问题 (Major)

#### 4. Fang 算法精度极差

- **位置**: `tdoa.c:tdoa_solve_fang()`
- **现象**: 误差133米（而Chan算法0.33米）
- **根因**: 
  - 默认值 `tij = 1.0` 当测量为0时（这本身就是问题）
  - TDOA测量获取逻辑有缺陷

```c
if (tij == 0.0) tij = 1.0;  // 硬编码默认值导致巨大误差
```

---

#### 5. Taylor 算法初始猜测不佳

- **位置**: `main.c:97`
- **现象**: Test 1 中 Taylor 误差2米，而 Chan 仅0.33米
- **根因**: 硬编码初始猜测 `{50.0, 28.0, 0.0}` 远离真实值 `{50.0, 30.0, 0.0}`

```c
tdoa_point3d_t initial_guess = {50.0, 28.0, 0.0};  // 应该用Chan结果
```

---

#### 6. 随机数种子问题

- **位置**: `main.c:130`
- **问题**: `srand(time(NULL))` 如果快速连续运行多次，种子相同

---

### 🟡 小问题 (Minor)

#### 7. Taylor 算法在 Test 5 中未收敛

- 误差20米，两个算法结果相同，可能陷入局部最优

#### 8. 缺少输入验证

- 未检查 `NaN` / `Inf` 结果
- 未检查迭代发散

#### 9. 打印格式问题

- `Test 3` 打印 `%.3f m` 但总时间是 `%.2f ms`，不一致

---

## 📈 代码质量评估

| 指标 | 评分 | 备注 |
|------|------|------|
| 编译警告 | 🟢 优秀 | 仅静态分析1个警告 |
| 内存管理 | 🟡 一般 | 存在潜在内存问题 |
| 错误处理 | 🟠 需改进 | 缺少边界检查 |
| 算法精度 | 🔴 需修复 | 多个算法有问题 |
| 代码注释 | 🟢 良好 | 文档较完整 |

---

## ✅ 测试断言建议

```c
// 1. 零噪声应该得到精确解
assert(err < 1e-6);

// 2. Monte Carlo 应该收敛到真实位置附近
assert(avg_error < 1.0);  // 1ns噪声对应约0.3米

// 3. 声学TDOA应该使用声速
assert(fabs(c - SPEED_OF_SOUND) < 1.0);  // 用于声学应用

// 4. 结果应该在合理范围内
assert(!isnan(result->x) && !isinf(result->x));
```

---

## 🎯 修复优先级

1. **立即修复**: Monte Carlo测试bug、声速常量问题
2. **高优先级**: Fang算法精度、初始猜测
3. **中优先级**: 内存安全问题、随机种子
4. **低优先级**: 打印格式、代码注释

---

**测试结论**: 代码基本可运行，但存在多个影响测试结果正确性的严重bug ⚠️
