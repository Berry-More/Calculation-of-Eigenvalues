#include <math.h>
#include <iostream>

long double PI = 3.14159265358979323846;

class Matrix
{

private:
	// matrix elements
	float
		x11, x12, x13,
		x21, x22, x23,
		x31, x32, x33;

	// matrix eigenvalues
	float
		l1, l2, l3;

public:

	// initialization function
	void init(float n11, float n12, float n13,
		float n21, float n22, float n23,
		float n31, float n32, float n33)
	{
		x11 = n11; x12 = n12; x13 = n13;
		x21 = n21; x22 = n22; x23 = n23;
		x31 = n31; x32 = n32; x33 = n33;
	}

	// print matrix in output
	void print() {
		std::cout << x11 << "\t" << x12 << "\t" << x13 << "\n";
		std::cout << x21 << "\t" << x22 << "\t" << x23 << "\n";
		std::cout << x31 << "\t" << x32 << "\t" << x33 << std::endl;
	}

	// matrix transpose
	void T() {
		init(x11, x21, x31, x12, x22, x32, x13, x23, x33);
	}

	// calculate determinant
	float getDeterminant() {
		return (x11 * (x22 * x33 - x23 * x32) 
			- x12 * (x21 * x33 - x23 * x31) 
			+ x13 * (x21 * x32 - x22 * x31));
	}

	// calculate eigenvalues
	// https://ru.intemodino.com/math/algebra/equations/cardano's-formula-for-solving-cubic-equations.html
	void eigvals() {

		// ax^3 + bx^2 + cx + d = 0

		float a, b, c, d;
		a = -1;
		b = x11 + x22 + x33;
		c = x12 * x21 + x13 * x31 - x11 * x33 
			- x11 * x22 - x22 * x33 + x23 * x32;
		d = x11 * x22 * x33 - x12 * x21 * x33 + x12 * x23 * x31
			+ x13 * x21 * x32 - x13 * x31 * x22 - x11 * x23 * x32;

		// I use Kardano equation
		// x = y - b/3a
		// y^3 + py + q = 0

		float p, q;
		p = (3 * a * c - pow(b, 2)) / (3 * pow(a, 2));
		q = (2 * pow(b, 3) - 9 * a * b * c + 27 * pow(a, 2) * d) 
			/ (27 * pow(a, 3));

		float q_big = pow((p / 3), 3) + pow((q / 2), 2);

		// if Q > 0: 1 real root and 2 complex roots
		if (q_big > 0) {

			float sign1 = (-q / 2 + pow(q_big, 0.5)) 
				/ abs(-q / 2 + pow(q_big, 0.5));
			float sign2 = (-q / 2 - pow(q_big, 0.5))
				/ abs(-q / 2 - pow(q_big, 0.5));

			float alpha = sign1 * pow(abs(-q / 2 + pow(q_big, 0.5)), 0.33);
			float betha = sign2 * pow(abs(-q / 2 - pow(q_big, 0.5)), 0.33);

			l1 = alpha + betha - b / (3 * a);
			l2 = (-(alpha + betha) / 2) - b / (3 * a);
			l3 = (alpha - betha) * pow(3, 0.5) / 2 ;

			std::cout << "l1 = " << l1 << std::endl;
			std::cout << "l2 = " << l2 << " + i(" << l3 << ")" << std::endl;
			std::cout << "l2 = " << l2 << " - i(" << l3 << ")" << std::endl;

		}

		// if Q = 0: 2 real roots
		if (q_big == 0) {

			l1 = 2 * pow((-q / 2), 0.33) - b / (3 * a);
			l2 = -1 * pow((-q / 2.), 0.33) - b / (3 * a);
			l3 = 0;

			std::cout << "l1 = " << l1 << std::endl;
			std::cout << "l2 = " << l2 << std::endl;
			std::cout << "l2 = " << l3 << std::endl;

		}

		// if Q < 0: 3 real roots
		if (q_big < 0) {

			float phi;

			if (q < 0) {
				phi = atan(pow(q_big, 0.5) / (-q / 2));
			}
			if (q > 0) {
				phi = atan(pow(q_big, 0.5) / (-q / 2)) + PI;
			}
			if (q == 0) {
				phi = PI / 2;
			}

			l1 = 2 * pow((p / 3), 0.5) * cos(phi / 3) - b / (3 * a);
			l2 = 2 * pow((p / 3), 0.5) * cos((phi / 3) + 2 * PI / 3) 
				- b / (3 * a);
			l3 = 2 * pow((p / 3), 0.5) * cos((phi / 3) + 4 * PI / 3) 
				- b / (3 * a);

			std::cout << "l1 = " << l1 << std::endl;
			std::cout << "l2 = " << l2 << std::endl;
			std::cout << "l2 = " << l3 << std::endl;

		}

	}

};

int main() {

	Matrix M1;

	M1.init(1, 0, 0, 0, 1, 0, 0, 0, 1);
	M1.print();

	std::cout << "\n" << std::endl;
	M1.eigvals();

	return 0;
}