#### 大整数

``` c++
#include<iostream>
#include<vector>
#include<deque>
#include<string.h>
#include<math.h>
#include<stdio.h>
#include<algorithm>
#include<string>
using namespace std;

class DividedByZeroException {};

class BigInteger {
private:
	vector<char> digits;
	bool sign;          //  true for positive, false for negitive
	void trim();        //  remove zeros in tail, but if the value is 0, keep only one:)
public:
	BigInteger(int);    // construct with a int integer
	BigInteger(string&);
	BigInteger();
	BigInteger(const BigInteger&);
	BigInteger operator=(const BigInteger& op2);

	BigInteger    abs() const;
	BigInteger    pow(int a);

	//binary operators

	friend BigInteger operator+=(BigInteger&, const BigInteger&);
	friend BigInteger operator-=(BigInteger&, const BigInteger&);
	friend BigInteger operator*=(BigInteger&, const BigInteger&);
	friend BigInteger operator/=(BigInteger&, const BigInteger&) throw(DividedByZeroException);
	friend BigInteger operator%=(BigInteger&, const BigInteger&) throw(DividedByZeroException);

	friend BigInteger operator+(const BigInteger&, const BigInteger&);
	friend BigInteger operator-(const BigInteger&, const BigInteger&);
	friend BigInteger operator*(const BigInteger&, const BigInteger&);
	friend BigInteger operator/(const BigInteger&, const BigInteger&) throw(DividedByZeroException);
	friend BigInteger operator%(const BigInteger&, const BigInteger&) throw(DividedByZeroException);


	//uniary operators
	friend BigInteger operator-(const BigInteger&);   //negative

	friend BigInteger operator++(BigInteger&);        //++v
	friend BigInteger operator++(BigInteger&, int);   //v++
	friend BigInteger operator--(BigInteger&);        //--v
	friend BigInteger operator--(BigInteger&, int);   //v--

	friend bool operator>(const BigInteger&, const BigInteger&);
	friend bool operator<(const BigInteger&, const BigInteger&);
	friend bool operator==(const BigInteger&, const BigInteger&);
	friend bool operator!=(const BigInteger&, const BigInteger&);
	friend bool operator>=(const BigInteger&, const BigInteger&);
	friend bool operator<=(const BigInteger&, const BigInteger&);

	friend ostream& operator<<(ostream&, const BigInteger&);   //print the BigInteger
	friend istream& operator>>(istream&, BigInteger&);         // input the BigInteger

public:
	static const BigInteger ZERO;
	static const BigInteger ONE;
	static const BigInteger TEN;
};
const BigInteger BigInteger::ZERO = BigInteger(0);
const BigInteger BigInteger::ONE = BigInteger(1);
const BigInteger BigInteger::TEN = BigInteger(10);


BigInteger::BigInteger() {
	sign = true;
}


BigInteger::BigInteger(int val) { // construct with a int integer
	if (val >= 0) {
		sign = true;
	}

	else {
		sign = false;
		val *= (-1);
	}

	do {
		digits.push_back((char)(val % 10));
		val /= 10;
	} while (val != 0);
}


BigInteger::BigInteger(string& def) {
	sign = true;

	for (string::reverse_iterator iter = def.rbegin(); iter < def.rend(); iter++) {
		char ch = (*iter);

		if (iter == def.rend() - 1) {
			if (ch == '+') {
				break;
			}

			if (ch == '-') {
				sign = false;
				break;
			}
		}

		digits.push_back((char)((*iter) - '0'));
	}

	trim();
}

void BigInteger::trim() {
	vector<char>::reverse_iterator iter = digits.rbegin();

	while (!digits.empty() && (*iter) == 0) {
		digits.pop_back();
		iter = digits.rbegin();
	}

	if (digits.size() == 0) {
		sign = true;
		digits.push_back(0);
	}
}


BigInteger::BigInteger(const BigInteger& op2) {
	sign = op2.sign;
	digits = op2.digits;
}


BigInteger BigInteger::operator=(const BigInteger& op2) {
	digits = op2.digits;
	sign = op2.sign;
	return (*this);
}


BigInteger BigInteger::abs() const {
	if (sign) {
		return *this;
	}

	else {
		return -(*this);
	}
}

BigInteger BigInteger::pow(int a) {
	BigInteger res(1);

	for (int i = 0; i < a; i++) {
		res *= (*this);
	}

	return res;
}

//binary operators
BigInteger operator+=(BigInteger& op1, const BigInteger& op2) {
	if (op1.sign == op2.sign) {     //只处理相同的符号的情况，异号的情况给-处理
		vector<char>::iterator iter1;
		vector<char>::const_iterator iter2;
		iter1 = op1.digits.begin();
		iter2 = op2.digits.begin();
		char to_add = 0;        //进位

		while (iter1 != op1.digits.end() && iter2 != op2.digits.end()) {
			(*iter1) = (*iter1) + (*iter2) + to_add;
			to_add = ((*iter1) > 9);    // 大于9进一位
			(*iter1) = (*iter1) % 10;
			iter1++;
			iter2++;
		}

		while (iter1 != op1.digits.end()) {    //
			(*iter1) = (*iter1) + to_add;
			to_add = ((*iter1) > 9);
			(*iter1) %= 10;
			iter1++;
		}

		while (iter2 != op2.digits.end()) {
			char val = (*iter2) + to_add;
			to_add = (val > 9);
			val %= 10;
			op1.digits.push_back(val);
			iter2++;
		}

		if (to_add != 0) {
			op1.digits.push_back(to_add);
		}

		return op1;
	}

	else {
		if (op1.sign) {
			return op1 -= (-op2);
		}

		else {
			return op1 = op2 - (-op1);
		}
	}

}

BigInteger operator-=(BigInteger& op1, const BigInteger& op2) {
	if (op1.sign == op2.sign) {     //只处理相同的符号的情况，异号的情况给+处理
		if (op1.sign) {
			if (op1 < op2) { // 2 - 3
				return  op1 = -(op2 - op1);
			}
		}

		else {
			if (-op1 > -op2) { // (-3)-(-2) = -(3 - 2)
				return op1 = -((-op1) - (-op2));
			}

			else {           // (-2)-(-3) = 3 - 2
				return op1 = (-op2) - (-op1);
			}
		}

		vector<char>::iterator iter1;
		vector<char>::const_iterator iter2;
		iter1 = op1.digits.begin();
		iter2 = op2.digits.begin();

		char to_substract = 0;  //借位

		while (iter1 != op1.digits.end() && iter2 != op2.digits.end()) {
			(*iter1) = (*iter1) - (*iter2) - to_substract;
			to_substract = 0;

			if ((*iter1) < 0) {
				to_substract = 1;
				(*iter1) += 10;
			}

			iter1++;
			iter2++;
		}

		while (iter1 != op1.digits.end()) {
			(*iter1) = (*iter1) - to_substract;
			to_substract = 0;

			if ((*iter1) < 0) {
				to_substract = 1;
				(*iter1) += 10;
			}

			else {
				break;
			}

			iter1++;
		}

		op1.trim();
		return op1;
	}

	else {
		if (op1 > BigInteger::ZERO) {
			return op1 += (-op2);
		}

		else {
			return op1 = -(op2 + (-op1));
		}
	}
}
BigInteger operator*=(BigInteger& op1, const BigInteger& op2) {
	BigInteger result(0);

	if (op1 == BigInteger::ZERO || op2 == BigInteger::ZERO) {
		result = BigInteger::ZERO;
	}

	else {
		vector<char>::const_iterator iter2 = op2.digits.begin();

		while (iter2 != op2.digits.end()) {
			if (*iter2 != 0) {
				deque<char> temp(op1.digits.begin(), op1.digits.end());
				char to_add = 0;
				deque<char>::iterator iter1 = temp.begin();

				while (iter1 != temp.end()) {
					(*iter1) *= (*iter2);
					(*iter1) += to_add;
					to_add = (*iter1) / 10;
					(*iter1) %= 10;
					iter1++;
				}

				if (to_add != 0) {
					temp.push_back(to_add);
				}

				int num_of_zeros = iter2 - op2.digits.begin();

				while (num_of_zeros--) {
					temp.push_front(0);
				}

				BigInteger temp2;
				temp2.digits.insert(temp2.digits.end(), temp.begin(), temp.end());
				temp2.trim();
				result = result + temp2;
			}

			iter2++;
		}

		result.sign = ((op1.sign && op2.sign) || (!op1.sign && !op2.sign));
	}

	op1 = result;
	return op1;
}

BigInteger operator/=(BigInteger& op1, const BigInteger& op2) throw(DividedByZeroException) {
	if (op2 == BigInteger::ZERO) {
		throw DividedByZeroException();
	}

	BigInteger t1 = op1.abs(), t2 = op2.abs();

	if (t1 < t2) {
		op1 = BigInteger::ZERO;
		return op1;
	}

	//现在 t1 > t2 > 0
	//只需将 t1/t2的结果交给result就可以了
	deque<char> temp;
	vector<char>::reverse_iterator iter = t1.digits.rbegin();

	BigInteger temp2(0);

	while (iter != t1.digits.rend()) {
		temp2 = temp2 * BigInteger::TEN + BigInteger((int)(*iter));
		char s = 0;

		while (temp2 >= t2) {
			temp2 = temp2 - t2;
			s = s + 1;
		}

		temp.push_front(s);
		iter++;
	}

	op1.digits.clear();
	op1.digits.insert(op1.digits.end(), temp.begin(), temp.end());
	op1.trim();
	op1.sign = ((op1.sign && op2.sign) || (!op1.sign && !op2.sign));
	return op1;
}

BigInteger operator%=(BigInteger& op1, const BigInteger& op2) throw(DividedByZeroException) {
	return op1 -= ((op1 / op2) * op2);
}

BigInteger operator+(const BigInteger& op1, const BigInteger& op2) {
	BigInteger temp(op1);
	temp += op2;
	return temp;
}
BigInteger operator-(const BigInteger& op1, const BigInteger& op2) {
	BigInteger temp(op1);
	temp -= op2;
	return temp;
}

BigInteger operator*(const BigInteger& op1, const BigInteger& op2) {
	BigInteger temp(op1);
	temp *= op2;
	return temp;

}

BigInteger operator/(const BigInteger& op1, const BigInteger& op2) throw(DividedByZeroException) {
	BigInteger temp(op1);
	temp /= op2;
	return temp;
}

BigInteger operator%(const BigInteger& op1, const BigInteger& op2) throw(DividedByZeroException) {
	BigInteger temp(op1);
	temp %= op2;
	return temp;
}

//uniary operators
BigInteger operator-(const BigInteger& op) {  //negative
	BigInteger temp = BigInteger(op);
	temp.sign = !temp.sign;
	return temp;
}

BigInteger operator++(BigInteger& op) {   //++v
	op += BigInteger::ONE;
	return op;
}

BigInteger operator++(BigInteger& op, int x) { //v++
	BigInteger temp(op);
	++op;
	return temp;
}

BigInteger operator--(BigInteger& op) {   //--v
	op -= BigInteger::ONE;
	return op;
}

BigInteger operator--(BigInteger& op, int x) { //v--
	BigInteger temp(op);
	--op;
	return temp;
}

bool operator<(const BigInteger& op1, const BigInteger& op2) {
	if (op1.sign != op2.sign) {
		return !op1.sign;
	}

	else {
		if (op1.digits.size() != op2.digits.size())
			return (op1.sign && op1.digits.size() < op2.digits.size())
			|| (!op1.sign && op1.digits.size() > op2.digits.size());

		vector<char>::const_reverse_iterator iter1, iter2;
		iter1 = op1.digits.rbegin();
		iter2 = op2.digits.rbegin();

		while (iter1 != op1.digits.rend()) {
			if (op1.sign &&  *iter1 < *iter2) {
				return true;
			}

			if (op1.sign &&  *iter1 > *iter2) {
				return false;
			}

			if (!op1.sign &&  *iter1 > *iter2) {
				return true;
			}

			if (!op1.sign &&  *iter1 < *iter2) {
				return false;
			}

			iter1++;
			iter2++;
		}

		return false;
	}
}
bool operator==(const BigInteger& op1, const BigInteger& op2) {
	if (op1.sign != op2.sign || op1.digits.size() != op2.digits.size()) {
		return false;
	}

	vector<char>::const_iterator iter1, iter2;
	iter1 = op1.digits.begin();
	iter2 = op2.digits.begin();

	while (iter1 != op1.digits.end()) {
		if (*iter1 != *iter2) {
			return false;
		}

		iter1++;
		iter2++;
	}

	return true;
}

bool operator!=(const BigInteger& op1, const BigInteger& op2) {
	return !(op1 == op2);
}

bool operator>=(const BigInteger& op1, const BigInteger& op2) {
	return (op1 > op2) || (op1 == op2);
}

bool operator<=(const BigInteger& op1, const BigInteger& op2) {
	return (op1 < op2) || (op1 == op2);
}

bool operator>(const BigInteger& op1, const BigInteger& op2) {
	return !(op1 <= op2);
}

ostream& operator<<(ostream& stream, const BigInteger& val) {  //print the BigInteger
	if (!val.sign) {
		stream << "-";
	}

	for (vector<char>::const_reverse_iterator iter = val.digits.rbegin(); iter != val.digits.rend(); iter++) {
		stream << (char)((*iter) + '0');
	}

	return stream;
}

istream& operator>>(istream& stream, BigInteger& val) {   //Input the BigInteger
	string str;
	stream >> str;
	val = BigInteger(str);
	return stream;
}

int main()
{
	string s = "0012032039243283298329";
	BigInteger a = s, b = 1000007;
	cout << a%b;
	return 0;
}
```
#### 素数测试

 - 前置： 快速乘、快速幂
 - int 范围内只需检查 2, 7, 61
 - long long 范围 2, 325, 9375, 28178, 450775, 9780504, 1795265022
 - 3E15内 2, 2570940, 880937, 610386380, 4130785767
 - 4E13内 2, 2570940, 211991001, 3749873356

``` c++
bool checkQ(LL a, LL n) {
    if (n == 2 || a >= n) return 1;
    if (n == 1 || !(n & 1)) return 0;
    LL d = n - 1;
    while (!(d & 1)) d >>= 1;
    LL t = bin(a, d, n);  // 不一定需要快速乘
    while (d != n - 1 && t != 1 && t != n - 1) {
        t = mul(t, t, n);
        d <<= 1;
    }
    return t == n - 1 || d & 1;
}

bool primeQ(LL n) {
    static vector <LL> t = {2, 325, 9375, 28178, 450775, 9780504, 1795265022};
    if (n <= 1) return false;
    for (LL k: t) if (!checkQ(k, n)) return false;
    return true;
}
```

#### 拓展欧几里得

``` c++
LL ex_gcd(LL a, LL b, LL &x, LL &y) {
    if (b == 0) { x = 1; y = 0; return a; }
    LL ret = ex_gcd(b, a % b, y, x);
    y -= a / b * x;
    return ret;
}
////////////////////////////////卡常版本
inline int ctz(LL x) { return __builtin_ctzll(x); }
LL gcd(LL a, LL b) {
    if (!a) return b; if (!b) return a;
    int t = ctz(a | b);
    a >>= ctz(a);
    do {
        b >>= ctz(b);
        if (a > b) swap(a, b);
        b -= a;
    } while (b);
    return a << t;
}
```

#### 逆元
- 如果p不是素数，使用拓展欧几里得
- 前置：快速幂/拓展欧几里得
``` c++
inline LL get_inv(LL x, LL p) { return bin(x, p - 2, p); }

LL get_inv(LL a, LL M) {
    static LL x, y;
    assert(exgcd(a, M, x, y) == 1);
    return (x % M + M) % M;
}

//• 预处理1~n 的逆元
LL inv[N];

void inv_init(LL n, LL p) {
    inv[1] = 1;
    FOR(i, 2, n)
    inv[i] = (p - p / i) * inv[p % i] % p;
}

//• 预处理阶乘及其逆元
LL invf[M], fac[M] = {1};

void fac_inv_init(LL n, LL p) {
    FOR(i, 1, n)
    fac[i] = i * fac[i - 1] % p;
    invf[n - 1] = bin(fac[n - 1], p - 2, p);
    FORD(i, n - 2, -1)
    invf[i] = invf[i + 1] * (i + 1) % p;
}
```

#### 组合数
 - 数较小，模数较大时使用逆元
 - 前置模板：逆元-预处理阶乘及其逆元

``` c++
inline LL c (LL n, LL m) { //n >= m >= 0
    return n < m || m < 0 ? 0 : fac[n] * invf[m] % MOD & invf[n-m] % mod
}
```
 - 预处理组合数
 ``` c++
 LL C[M][M];
void init_C(int n) {
    FOR (i, 0, n) {
        C[i][0] = C[i][i] = 1;
        FOR (j, 1, i)
        C[i][j] = (C[i - 1][j] + C[i - 1][j - 1]) % MOD;
    }
}
 ```
 
 #### 素数筛
 ```c++
 const int maxn = 1e5 + 50;
int prime[maxn];
bool is_prime[maxn];
int sieve(int n) {
	int p = 0;
	memset(is_prime, true, sizeof(is_prime));
	is_prime[0] = is_prime[1] = false;
	for (int i = 2; i <= n; i++) {
		if (is_prime[i]) {
			prime[p++] = i;
			for (int j = 2 * i; j <= n; j+=i) {
				is_prime[j] = false;
			}
		}
	}
	return p;
}
 ```
 #### 中国剩余定理
 ```c++
 LL CRT(LL *m, LL *r, LL n) {
	if (!n) return 0;
	LL M = m[0], R = r[0], x, y, d;
	FOR (i, 1, n) {
		d = ex_gcd(M, m[i], x, y);
		if ((r[i] - R) % d) return -1;
		x = (r[i] - R) / d * x % (m[i] / d);
		R += x * M;
		M = M / d * m[i];
		R %= M;
	}
	return R >= 0 ? R : R + M;
}
 ```
 
 #### 线性基
 ```c++
 /*
 * 多个线性基之间可以进行合并，即
 * 把一个线性基中的元素全部插入另一个中。
 */
struct Linear_Basis
{
    LL b[63], nb[63], tot;

    void init()
    {
        tot = 0;
        memset(b, 0, sizeof(b));
        memset(nb, 0, sizeof(nb));
    }

    bool ins(LL x)//插入
    {
        for (int i = 62; i >= 0; i--)
            if (x&(1LL << i))
            {
                if (!b[i]) { b[i] = x; break; }
                x ^= b[i];
            }
        return x > 0;
    }

    LL Max(LL x)//所有可能异或中最大
    {
        LL res = x;
        for (int i = 62; i >= 0; i--)
            res = max(res, res^b[i]);
        return res;
    }

    LL Min(LL x)//所有可能异或中最小
    {
        LL res = x;
        for (int i = 0; i <= 62; i++)
            if (b[i]) res ^= b[i];
        return res;
    }

    void rebuild()
    {
        for (int i = 62; i >= 0; i--)
            for (int j = i - 1; j >= 0; j--)
                if (b[i] & (1LL << j)) b[i] ^= b[j];
        for (int i = 0; i <= 62; i++)
            if (b[i]) nb[tot++] = b[i];
    }

    LL Kth_Max(LL k) //所有可能异或中k小
    {
        rebuild();
        LL res = 0;
        for (int i = 62; i >= 0; i--)
            if (k&(1LL << i)) res ^= nb[i];
        return res;
    }

} LB;
 ```
 
 #### 博弈
·Nim 游戏：每轮从若干堆石子中的一堆取走若干颗。先手必胜条件为石子数量异或和非零。
·阶梯 Nim 游戏：可以选择阶梯上某一堆中的若干颗向下推动一级，直到全部推下去。先手必胜条件是奇数阶梯的异或和非零（对于偶数阶梯的操作可以模仿）。
 - Anti-SG：无法操作者胜。先手必胜的条件是：
 - SG 不为 0 且某个单一游戏的 SG 大于 1 。
 - SG 为 0 且没有单一游戏的 SG 大于 1。
 - Every-SG：对所有单一游戏都要操作。先手必胜的条件是单一游戏中的最大 step 为奇数。
 - 对于终止状态 step 为 0
 - - 对于 SG 为 0 的状态，step 是最大后继 step +1
 - - 对于 SG 非 0 的状态，step 是最小后继 step +1
 - 树上删边：叶子 SG 为 0，非叶子结点为所有子结点的 SG 值加 1 后的异或和。
 - 尝试：
 - - 打表找规律
 - - 寻找一类必胜态（如对称局面）
 - - 直接博弈 dp
