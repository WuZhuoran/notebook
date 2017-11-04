import java.util.*;
import javafx.util.*;

public class MathKibrary {

	// Greatest Common Divisor
	public static int gcd(int a, int b) {
		return b == 0 ? a : gcd(b, a % b);
	}

	// Least Common Multiple
	public static int lcm(int a, int b) {
		return a * (b / gcd(a, b));
	}

	// Prime Numbers
	public static int _sieve_size;
	public static boolean[] bs; // 10^7 should be enough for most cases
	public static List<Integer> primes = new ArrayList<Integer>(); // compact
																	// list of
																	// primes

	// first part

	public static void sieve(int upperbound) { // create list of primes in
												// [0..upperbound]
		_sieve_size = upperbound + 1; // add 1 to include upperbound
		bs = new boolean[_sieve_size];
		Arrays.fill(bs, true); // set all bits to 1
		bs[0] = bs[1] = false; // except index 0 and 1
		for (long i = 2; i < _sieve_size; i++) {
			if (bs[(int) i]) {
				// cross out multiples of i starting from i * i!
				for (long j = i * i; j < _sieve_size; j += i) {
					bs[(int) j] = false;
				}
				primes.add((int) i); // also add this vector containing list of
										// primes
			}
		}
	} // call this method in main method

	public static boolean isPrime(long N) { // a good enough deterministic prime
											// tester
		if (N < _sieve_size) {
			return bs[(int) N]; // O(1) for small primes
		}
		for (int i = 0; i < primes.size(); i++) {
			if (N % primes.get(i) == 0) {
				return false;
			}
		}
		return true; // it takes longer time if N is a large prime!
	} // note: only work for N <= (last prime in vi "primes")^2

	// Find Prime Factors
	public static List<Integer> primeFactors(long N) {
		List<Integer> factors = new ArrayList<Integer>();
		int PF_idx = 0;
		long PF = primes.get(PF_idx);
		while (PF * PF <= N) {
			while (N % PF == 0) {
				N /= PF;
				factors.add((int) PF);
			}
			PF = primes.get(++PF_idx); // only consider primes
		}
		if (N != 1) {
			factors.add((int) N);
		}

		return factors;
	}

	// Count the number of prome factors of N
	public static long numPF(long N) {
		int PF_idx = 0;
		long PF = primes.get(PF_idx);
		long ans = 0L;
		while (PF * PF <= N) {
			while (N % PF == 0) {
				N /= PF;
				ans++;
			}
			PF = primes.get(++PF_idx);
		}
		if (N != 1) {
			ans++;
		}
		return ans;
	}

	// Cpimnt the number of different prime factors of N
	public static long numDiffPF(long N) {
		int PF_idx = 0;
		long PF = primes.get(PF_idx);
		long ans = 0;
		while (N != 1 && (PF * PF <= N)) {
			if (N % PF == 0) {
				ans++; // count this pf only once
			}
			while (N % PF == 0) {
				N /= PF;
			}
			PF = primes.get(++PF_idx);
		}
		if (N != 1) {
			ans++;
		}
		return ans;
	}

	// Sum the prime factors of N
	public static long sumPF(long N) {
		int PF_idx = 0;
		long PF = primes.get(PF_idx);
		long ans = 0;
		while (N != 1 && (PF * PF <= N)) {
			while (N % PF == 0) {
				N /= PF;
				ans += PF;
			}
			PF = primes.get(++PF_idx);
		}
		if (N != 1) {
			ans += N;
		}
		return ans;
	}

	// Count the number of divisors
	public static long numDiv(long N) {
		int PF_idx = 0;
		long PF = primes.get(PF_idx);
		long ans = 1;
		while (N != 1 && (PF * PF <= N)) {
			long power = 0;
			while (N % PF == 0) {
				N /= PF;
				power++;
			}
			ans *= (power + 1);
			PF = primes.get(++PF_idx);
		}

		if (N != 1) {
			ans *= 2;
		}
		return ans;
	}

	// Summ the divisors
	public static long sumDiv(long N) {
		int PF_idx = 0;
		long PF = primes.get(PF_idx);
		long ans = 1;
		while (N != 1 && (PF * PF <= N)) {
			long power = 0;
			while (N % PF == 0) {
				N /= PF;
				power++;
			}
			ans *= ((long) Math.pow((double) PF, power + 1.0) - 1) / (PF - 1);
			PF = primes.get(++PF_idx);
		}

		if (N != 1) {
			ans *= ((long) Math.pow((double) N, 2.0) - 1) / (N - 1);
		}
		return ans;
	}

	// EulerPhi: Count the number of positive integers < N that are relatively
	// prime to N
	public static long EulerPhi(long N) {
		int PF_idx = 0;
		long PF = primes.get(PF_idx);
		long ans = N;

		while (N != 1 && (PF * PF <= N)) {
			if (N % PF == 0) {
				ans -= ans / PF;
			}
			while (N % PF == 0) {
				N /= PF;
			}
			PF = primes.get(++PF_idx);
		}

		if (N != 1) {
			ans -= ans / N;
		}
		return ans;
	}

	public static int x;
	public static int y;
	public static int d;

	// store x, y and d as global variables
	public static void extendedEuclid(int a, int b) {
		if (b == 0) {
			x = 1;
			y = 0;
			d = a;
			return;
		}
		extendedEuclid(b, a % b);
		int x1 = y;
		int y1 = x - (a / b) * y;
		x = x1;
		y = y1;
	}

	public static int f(int x) {
		return (3 * x + 1) % 4;
	}

	public static Pair<Integer, Integer> floydCycleFinding(int x0) {
		// 1st part, finding k * mu, hare's speed is 2x tortoise's
		int tortoise = f(x0);
		int hare = f(f(x0));
		while (tortoise != hare) {
			tortoise = f(tortoise);
			hare = f(f(hare));
		}

		// 2nd part, finding mu, hare and tortoise move at the same speed
		int mu = 0;
		hare = x0;
		while (tortoise != hare) {
			tortoise = f(tortoise);
			hare = f(hare);
			mu++;
		}

		// 3rd part, finding lambda, hare moves, tortoise stays
		int lambda = 1;
		hare = f(tortoise);
		while (tortoise != hare) {
			hare = f(hare);
			lambda++;
		}

		return new Pair<Integer, Integer>(mu, lambda);
	}

	public static void main(String[] args) {
		sieve(100000000);
		System.out.println(isPrime(2147483647));
		System.out.println(isPrime(136117223861L));

		List<Integer> r = primeFactors(2147483647);
		for (int i = 0; i < r.size(); i++) {
			System.out.println(r.get(i));
		}
		r = primeFactors(136117223861L);
		for (int i = 0; i < r.size(); i++) {
			System.out.println(r.get(i));
		}

		System.out.printf("numPF(%d) = %d\n", 50, numPF(50)); // 2^1 * 5^2 => 3
		System.out.printf("numDiffPF(%d) = %d\n", 50, numDiffPF(50)); // 2^1 *
																		// 5^2
																		// => 2
		System.out.printf("sumPF(%d) = %d\n", 50, sumPF(50)); // 2^1 * 5^2 => 2
																// + 5 + 5 = 12
		System.out.printf("numDiv(%d) = %d\n", 50, numDiv(50)); // 1, 2, 5, 10,
																// 25, 50, 6
																// divisors
		System.out.printf("sumDiv(%d) = %d\n", 50, sumDiv(50)); // 1 + 2 + 5 +
																// 10 + 25 + 50
																// = 93
		System.out.printf("EulerPhi(%d) = %d\n", 50, EulerPhi(50)); // 20
																	// integers
																	// < 50 are
																	// relatively
																	// prime
																	// with 50

	}

}
