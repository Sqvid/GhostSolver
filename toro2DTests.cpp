double explLS(double x, double y) {
	double r = 0.4;

	return x*x + y*y - r*r;
}

double cylExplDensity(double x, double y) {
	double phi = explLS(x, y);

	return phi <= 0 ? 1 : 0.125;
}

double cylExplVelocityX(double x, double y) {
	x = 0;
	y = 0;

	return 0;
}

double cylExplVelocityY(double x, double y) {
	x = 0;
	y = 0;

	return 0;
}

double cylExplPressure(double x, double y) {
	double phi = explLS(x, y);

	return phi <= 0 ? 1 : 0.1;
}
