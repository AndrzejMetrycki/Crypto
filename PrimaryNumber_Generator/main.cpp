#include <stdint.h>
#include <iostream>

#include <thread>
#include <vector>
#include <mutex>
#include <chrono>
#include <cstdint>
#include <math.h>



using namespace std;
using namespace chrono;

//__int128 mint; 
typedef uint64_t _int_t;

class CPrimeNr_Virtual
{
public:
	virtual bool IsPrime(_int_t iNr) = 0;
};




class CPrimeNr : public CPrimeNr_Virtual
{
public:


	virtual void OnPrimeNumberFound(_int_t iNr)
	{
		static mutex mMutex;
		static int iCounter = 0;

		mMutex.lock();
		cout << "Liczba Pierwsza[" << iCounter++ << "]: " << iNr << endl;
		mMutex.unlock();

	}

	void FindPrimeNumbers(_int_t iBegin, _int_t iEnd)
	{
		for (_int_t i = iBegin; i <= iEnd; ++i)
		{
			if (IsPrime(i))
			{
				OnPrimeNumberFound(i);
			}
		}
	}



	void Multithread_FindPrimaryNumbers(_int_t iBegin, _int_t iEnd)
	{
		vector<thread> vThd;
		uint64_t iTotalThdCount = 32;
		if ((iEnd - iBegin) < iTotalThdCount)
			iTotalThdCount = iEnd - iBegin;

		uint64_t iCounter = (iEnd - iBegin) / iTotalThdCount;



		for (uint64_t i = 0; i != iTotalThdCount; ++i)
		{
			uint64_t iNext = iBegin + iCounter;
			thread thd{ &CPrimeNr::FindPrimeNumbers,this ,iBegin, iNext };
			vThd.push_back(move(thd));
			iBegin = iNext;
		}
		if (iBegin != iEnd)
		{
			thread thd{ &CPrimeNr::FindPrimeNumbers,this ,iBegin, iEnd };
			vThd.push_back(move(thd));
		}

		for (uint32_t i = 0; i != vThd.size(); ++i)
		{
			vThd[i].join();
		}

	}
};



class CPrimeNr_SimpleMethod : public CPrimeNr
{
public:
	bool IsPrime(_int_t iNr)
	{
		if (iNr <= 1)
			return false;

		if (iNr <= 3)
			return true;
		else if (!(iNr % 2) || !(iNr % 3))
			return false;

		

		_int_t iEnd = (_int_t)sqrt(iNr) + 1;

		for (_int_t i = 5;i <= iEnd;i += 6)
		{
			if (!(iNr%i) || !(iNr % (i + 2)))
				return false;
		}
		return true;
	}
};


class CPrimeNr_Fermat : public CPrimeNr
{
public:
	CPrimeNr_Fermat(void)
	{
		srand(time(0));
	}
	bool IsPrime(_int_t iNr)
	{
		return Fermat(iNr, 50);
	}

private:

	_int_t modulo(_int_t base, _int_t exponent, _int_t mod)
	{
		_int_t x = 1;
		_int_t y = base;
		while (exponent > 0)
		{
			if (exponent % 2 == 1)
				x = (x * y) % mod;
			y = (y * y) % mod;
			exponent = exponent / 2;
		}
		return x % mod;
	}

	bool Fermat(_int_t p, int iterations)
	{
		if (p == 1)
		{
			return false;
		}
		for (int i = 0; i < iterations; i++)
		{
			_int_t a = rand() % (p - 1) + 1;
			if (modulo(a, p - 1, p) != 1)
			{
				return false;
			}
		}
		return true;
	}
};

class CPrimeNr_Miller : public CPrimeNr
{
public:
	CPrimeNr_Miller(void)
	{
		srand(time(0));
	}
	
	_int_t modulo(_int_t base, _int_t exponent, _int_t mod)
	{
		_int_t x = 1;
		_int_t y = base;
		while (exponent > 0)
		{
			if (exponent % 2 == 1)
				x = (x * y) % mod;
			y = (y * y) % mod;
			exponent = exponent / 2;
		}
		return x % mod;
	}
	_int_t mulmod(_int_t a, _int_t b, _int_t mod)
	{
		_int_t x = 0, y = a % mod;
		while (b > 0)
		{
			if (b % 2 == 1)
			{
				x = (x + y) % mod;
			}
			y = (y * 2) % mod;
			b /= 2;
		}
		return x % mod;
	}

	bool IsPrime(_int_t n)
	{
		if (n <= 1)
			return false;

		if (n <= 3)
			return true;

		if (!(n % 2))
			return false;

		register _int_t d = n - 1;
		register _int_t r = 0;

		while (!(d&1))
		{
			d >>= 1;
			r++;
		}


		register int k = 50;

		for (register int i = 0; i != k; ++i)
		{
			_int_t a = rand() % (n - 3) + 2;
			//_int_t a = rand() % (n - 1) + 1;
			//_int_t x = Pow(a, d) % n;
			_int_t x = modulo(a, d, n);
/*
			if (x1 != x)
			{
				int dupa = 55;
				dupa++;
			}*/

			if ((x == 1) || (x == (n-1)))
				continue;

			if (InLoop(x, n,r))
				continue;

			return false;
		}
		return true;
	}

private:
	bool InLoop(_int_t &x, _int_t &n, register _int_t r)
	{
		while (--r)
		{
			_int_t x1 = mulmod(x, x, n);

			x = (x*x) % n;

			if (x1 != x)
			{
				int dupa = 55;
				dupa++;
			}


			if (x == 1)
				return false;

			if (x == n-1)
				return true;
		}
		return false;
	}
	_int_t Pow(_int_t x, _int_t n)
	{

		if (!n)
			return 1;

		if (n&1)
		{
			_int_t iTmp = Pow(x, ((n - 1)/2));
			return x*iTmp*iTmp;
		}
		else
		{
			_int_t iTmp = Pow(x, (n / 2));
			return iTmp*iTmp;
		}
	}


};

class CPrimeNr_Rabin : public CPrimeNr
{
public:
	CPrimeNr_Rabin(void)
	{
		srand(time(0));
	}

	using ll = _int_t;

	bool IsPrime(_int_t iNr)
	{
		return Miller(iNr, 50);
	}


private:
	ll mulmod(ll a, ll b, ll mod)
	{
		ll x = 0, y = a % mod;
		while (b > 0)
		{
			if (b % 2 == 1)
			{
				x = (x + y) % mod;
			}
			y = (y * 2) % mod;
			b /= 2;
		}
		return x % mod;
	}
	/*
	 * modular exponentiation
	  */
	ll modulo(ll base, ll exponent, ll mod)
	{
		ll x = 1;
		ll y = base;
		while (exponent > 0)
		{
			if (exponent % 2 == 1)
				x = (x * y) % mod;
			y = (y * y) % mod;
			exponent = exponent / 2;
		}
		return x % mod;
	}

	/*
	 * Miller-Rabin primality test, iteration signifies the accuracy
	  */
	bool Miller(ll p, int iteration)
	{
		if (p < 2)
		{
			return false;
		}
		if (p != 2 && p % 2 == 0)
		{
			return false;
		}
		ll s = p - 1;
		while (s % 2 == 0)
		{
			s /= 2;
		}
		for (int i = 0; i < iteration; i++)
		{
			ll a = rand() % (p - 1) + 1, temp = s;
			ll mod = modulo(a, temp, p);
			while (temp != p - 1 && mod != 1 && mod != p - 1)
			{
				mod = mulmod(mod, mod, p);
				temp *= 2;
			}
			if (mod != p - 1 && temp % 2 == 0)
			{
				return false;
			}
		}
		return true;
	}
};


bool IsPrimaryNumber(uint64_t iNr)
{
	if (iNr < 2)
		throw string("Wrong Nr");

	if (!(iNr % 2))
		return false;

	uint64_t iEnd = sqrt((double)iNr)+1;
	for (uint64_t i = 3; i <= iEnd; i+=2)
	{
		if (!(iNr%i))
			return false;
	}
	return true;
}

void OnPrimaryNumberFound(uint64_t iNr)
{
	static mutex mMutex;
	static int iCounter = 0;

	mMutex.lock();
	cout << "Liczba Pierwsza[" << iCounter++ <<"]: "<< iNr << endl;
	mMutex.unlock();

}

void FindPrimaryNumbers(uint64_t iBegin, uint64_t iEnd)
{
	if (iBegin < 2)
		iBegin = 3;
	if (iBegin > iEnd)
		return;

	if (!(iBegin % 2))
		iBegin++;

	for (uint64_t i = iBegin;i <= iEnd; i+=2)
	{
		if (IsPrimaryNumber(i))
		{
			OnPrimaryNumberFound(i);
		}
	}
}


void MultiThread_FindPrimaryNumbers(uint64_t iBegin, uint64_t iEnd)
{
	vector<thread> vThd;
	uint64_t iTotalThdCount = 16;
	if ((iEnd - iBegin) < iTotalThdCount)
		iTotalThdCount = iEnd - iBegin;

	uint64_t iCounter = (iEnd - iBegin) / iTotalThdCount;
	for (uint64_t i = 0;i != iTotalThdCount; ++i)
	{
		uint64_t iNext = iBegin + iCounter;
		thread thd{ FindPrimaryNumbers,iBegin, iNext };
		vThd.push_back(move(thd));
		iBegin = iNext;
	}
	if (iBegin != iEnd)
	{
		thread thd{ FindPrimaryNumbers,iBegin, iEnd };
		vThd.push_back(move(thd));
	}

	for (uint32_t i = 0; i != vThd.size(); ++i)
	{
		vThd[i].join();
	}


}


class CPrime_Strong
{
public:
	bool IsPrime(_int_t iNr)
	{
		CPrimeNr_SimpleMethod mMth;

		if (mMth.IsPrime(iNr))
		{
			iNr = (iNr-1)/2;
			
			return mMth.IsPrime(iNr);

		}
		return false;
	}
};

class CFactoring_Simple
{
public:

	vector<_int_t> Factor(_int_t iNr)
	{
		vector<_int_t> viTmp;
		_int_t iResp = FactorBase(iNr);
		viTmp.push_back(iResp);
		if (iResp == iNr)
		{
			return viTmp;
		}
		else
		{
			auto viResp = Factor(iNr / iResp);

			for (int i = 0; i != viResp.size(); ++i)
			{
				viTmp.push_back(viResp[i]);
			}
		}

		return viTmp;
	}




private:
	_int_t FactorBase(_int_t iNr)
	{
		_int_t iEnd = (_int_t)sqrt(iNr) + 1;
		while (iEnd*iEnd < iNr)
			iEnd++;

		for (_int_t i = 5; i <= iEnd; i += 6)
		{
			if (!(iNr%i))
				return i;
				
			if (!(iNr % (i + 2)))
			{
				return i + 2;
				break;
			}
		}

		return iNr;
	}
};

void Factoring(void)
{
	CPrime_Strong mSimple;
	vector<_int_t> m_viPrime;
	for (_int_t i = 100000000; ; ++i)
	{
		if (mSimple.IsPrime(i))
		{
			m_viPrime.push_back(i);
			if (m_viPrime.size() == 2)
				break;
		}
	}

	cout << "Wygenerowane liczby" << endl;
	for (uint32_t i = 0; i != m_viPrime.size(); ++i)
	{
		cout << "L[" << i << "]: " << m_viPrime[i] << endl;
	}
	_int_t iMulNr = m_viPrime[0] * m_viPrime[1];
	//iMulNr += 2;
	cout << "Iloczyn liczb: " << iMulNr <<endl;

	CFactoring_Simple mFactor;

	auto t1 = steady_clock::now();
	auto viFactor = mFactor.Factor(iMulNr);


	
	auto d = steady_clock::now() - t1; 
	auto ms = duration_cast<milliseconds>(d).count();

	cout << "Czas: " <<(double)ms/1000 <<endl;


	for (uint32_t i = 0; i != viFactor.size(); ++i)
	{
		cout << "LF[" << i << "]: " << viFactor[i] << endl;
	}



}



int main(void)
{
	{
		CD tmp;
	}



	Factoring();
	system("pause");
	return 0;


	cout << "Czesc" << endl;
	uint64_t iBegin = 2305843009213693951;
	uint64_t iEnd = iBegin + 160;

	CPrimeNr_SimpleMethod mSimple;
	CPrimeNr_Miller mMiller;
	CPrimeNr_Rabin mRabin;

	


	for (_int_t i = 2; i != 100000; ++i)
	{
		bool b1 = mSimple.IsPrime(i);
		bool b2 = mMiller.IsPrime(i);
		bool b3 = mRabin.IsPrime(i);
		if (b1 != b2 || b3 != b2)
		{
			cout << "Rozbierznoœæ: " << i << "\t\t B1:" << b1 << " \tB2: " << b2 << " \tB3: " << b3 << endl;
		}


	}
	system("pause");
	return 0;

	/*
	for (int x = 1; x != 20; ++x)
	{
		for (int y = 0; y != 10; ++y)
		{
			_int_t iRes1 = Pow(x,y);
			_int_t iRes2 = pow(x, y);

			cout <<"x:"<<x<<"y: "<<y <<"\t\t" << "Res1: " << iRes1 << "\t\t" << "Res2: " << iRes2 << endl;

			if (iRes1 != iRes2)
			{
				system("pause");
				return 0;
			}

		}
	}
	system("pause");
	return 0;
*/
	cout << "Fast search" << endl;
	
	steady_clock::time_point t1;
	int64_t ms;
	steady_clock::duration d;
	t1 = steady_clock::now();
	MultiThread_FindPrimaryNumbers(iBegin, iEnd);
	d = steady_clock::now() - t1; 	ms = duration_cast<milliseconds>(d).count();
	cout << "TimeSpend[s]: " << (double)ms / 1000.0f << endl;
/*
	cout << "Slow search" << endl;
	t1 = steady_clock::now();
	FindPrimaryNumbers(iBegin, iEnd);
	d = steady_clock::now() - t1; 	ms = duration_cast<milliseconds>(d).count();
	cout << "TimeSpend[s]: " << (double)ms / 1000.0f << endl;
	*/

/*
	cout << "Class search" << endl;
	t1 = steady_clock::now();
	CPrimeNr mPrimeNr;
	mPrimeNr.Multithread_FindPrimaryNumbers(iBegin, iEnd);
	d = steady_clock::now() - t1; 	ms = duration_cast<milliseconds>(d).count();
	cout << "TimeSpend[s]: " << (double)ms / 1000.0f << endl;

	*/
	system("pause");
	return 0;
}
