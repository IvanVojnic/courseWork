package main

import (
	"fmt"
	"github.com/sbinet/go-gnuplot"
	"math"
)

func infinityNorm(a []float64, b []float64) float64 {
	if len(a) != len(b) {
		panic("Vectors must be of the same length.")
	}

	maxDiff := 0.0
	for i := 0; i < len(a); i++ {
		diff := math.Abs(a[i] - b[i])
		if diff > maxDiff {
			maxDiff = diff
		}
	}
	return maxDiff
}

func kapSecond() {

	var rList [][]float64
	var zList [][]float64

	Rho := 1.0                  // плотность жидкости, г/см^3
	Sigma := 72.75              // коэффициент поверхностного натяжения, дин/см
	V := 1.0                    // объем капли, см^3
	angle := 45 * math.Pi / 180 // угол смачивания, градусы
	n := 1000                   // количество разбиений
	h := 1.0 / float64(n)       // величина шага
	eps := math.Pow(h, 5)       // допустимая погрешность
	g := 981.0                  // ускорение свободного падения, см/сек^2
	Bo := 0.0                   // Rho*g/(Sigma*math.Pow(V, 0.66))
	z := make([]float64, n+1)
	zprev := make([]float64, n+1)
	r := make([]float64, n+1)
	a := make([]float64, n+1)
	A := make([]float64, n-1)
	B := make([]float64, n-1)
	C := make([]float64, n-1)
	F := make([]float64, n-1)
	alpha := make([]float64, n)
	beta := make([]float64, n)

	zRazm := make([]float64, n+1)
	rRazm := make([]float64, n+1)

	for i := 0; i < n-1; i++ {
		A[i] = 0
		B[i] = 0
		C[i] = 0
		F[i] = 0
	}

	// начальное приближение
	for i := 0; i < n+1; i++ {
		r[i] = float64(i) * h
		z[i] = 1 - r[i]
	}

	q := 0
	sum := 0.0
	tau := 0.5
	I := 0.0
	Q := 0.0

	m1 := 0.0
	v1 := 0.0
	v2 := 0.0

	for k := 0; k < 3; k++ {
		switch k {
		case 0:
			Bo = 0
		case 1:
			Bo = 1
		case 2:
			Bo = 2
		}

		zprev = make([]float64, n+1)
		q = 0

		for infinityNorm(zprev, z) > (1/tau)*eps && q < 1000 {

			copy(zprev, z)

			a[0] = 0
			for i := 1; i <= n-1; i++ {
				a[i] = 0.5 * (r[i-1] + r[i]) / math.Sqrt(1+math.Pow((z[i]-z[i-1])/h, 2))
			}
			a[n] = 0

			sum = 0
			for i := 1; i <= n-1; i++ {
				sum += r[i] * z[i]
			}

			I = 2 * math.Pi * h * sum
			Q = -2*math.Sin(angle) - (Bo * math.Pow(I, 1.0/3.0) / math.Pi)

			for i := 1; i <= n-2; i++ {
				A[i] = a[i]
			}

			for i := 1; i <= n-2; i++ {
				B[i] = a[i+1]
			}

			for i := 1; i <= n-2; i++ {
				C[i] = a[i+1] + a[i] + Bo*r[i]*math.Pow(h, 2)/math.Pow(I, 2.0/3.0)
			}

			for i := 1; i <= n-2; i++ {
				F[i] = -r[i] * math.Pow(h, 2) * Q
			}

			m1 = 1 / (1 + (math.Pow(h, 2)*Bo)/4)
			v1 = -(math.Pow(h, 2) * Q) / (4 * (1 + (math.Pow(h, 2)*Bo)/4))

			v2 = h*math.Tan(angle) + (math.Pow(h, 2) * (Q + math.Sin(angle)) / (2 * math.Pow(math.Cos(angle), 3)))

			alpha[1] = m1
			beta[1] = v1

			for i := 1; i <= n-2; i++ {
				alpha[i+1] = B[i] / (C[i] - (alpha[i] * A[i]))
				beta[i+1] = ((A[i] * beta[i]) + F[i]) / (C[i] - (alpha[i] * A[i]))
			}

			z[n] = 0
			z[n-1] = v2

			for i := n - 2; i >= 0; i-- {
				z[i] = alpha[i+1]*z[i+1] + beta[i+1]
			}

			for i := 0; i < n+1; i++ {
				z[i] = tau*z[i] + (1-tau)*zprev[i]
			}
			q++
		}

		if q >= 1000 {
			fmt.Println("Не удалось достичь заданной точности за 1000 итераций:", q)
		} else {
			fmt.Println("Успешно достигнута заданная точность за", q, "итераций")
		}

		fmt.Println("Параметры расчета:")
		fmt.Println("Плотность жидкости:", Rho, "г/см^3")
		fmt.Println("Коэффициент поверхностного натяжения:", Sigma, "дин/см")
		fmt.Println("Ускорение свободого падения:", g, "см/сек^2")
		fmt.Println("Объем капли:", V, "см^3")
		fmt.Println("Угол смачивания:", angle*180/math.Pi, "градусов")
		fmt.Println("Коэффициент Бонджера (Bo):", Bo)

		// пересчет координат точек для размещения на графике
		for i := 0; i < n+1; i++ {
			zRazm[i] = z[i] * 1000
			rRazm[i] = r[i] * 1000
		}

		// добавление значений в списки для построения графика
		rList = append(rList, rRazm)
		zList = append(zList, zRazm)
	}

	// построение графика
	plt, _ := gnuplot.NewPlotter("", false, false)
	defer plt.Close()

	plt.CheckedCmd("set term wxt size 700,500")
	plt.CheckedCmd("set title 'Капиллярный контур на поверхности жидкости'")
	plt.CheckedCmd("set xlabel 'Радиус, мкм'")
	plt.CheckedCmd("set ylabel 'Высоть, мкм'")

	plt.CheckedCmd("set key top left")
	plt.CheckedCmd("set grid")

	plt.PlotXY(rList[0], zList[0], "Коэффициент Бонджера (Bo) = 0")
	plt.PlotXY(rList[1], zList[1], "Коэффициент Бонджера (Bo) = 10")
	plt.PlotXY(rList[2], zList[2], "Коэффициент Бонджера (Bo) = 100")
	plt.CheckedCmd("set xrange [-10:510]")
	plt.CheckedCmd("set yrange [-10:510]")

	plt.CheckedCmd("pause -1")
}

func main() {
	kapSecond()
}
