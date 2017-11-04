import java.util.*;

public class GeometryLibrary {
	public static final double EPS = 1e-9;

	public static final double PI = Math.PI;

	public static final double PI_ACOS = Math.acos(-1.0);

	/**
	 * Points
	 */

	public static class Point_i implements Comparable<Point_i> {
		public int x;
		public int y;

		Point_i() {
			x = 0;
			y = 0;
		}

		Point_i(int x, int y) {
			this.x = x;
			this.y = y;
		}

		@Override
		public int compareTo(Point_i o) {
			// Ascending
			return (x == o.x) ? y - o.y : x - o.x;
			// Descending
			// return (x == o.x) ? o.y - y : o.x - x;
		}

		public boolean equals(Point_i obj) {
			return (x == obj.x) && (y == obj.y);
		}
	}

	public static class Point_d implements Comparable<Point_d> {
		public double x;
		public double y;

		Point_d() {
			x = 0.0;
			y = 0.0;
		}

		Point_d(double x, double y) {
			this.x = x;
			this.y = y;
		}

		@Override
		public int compareTo(Point_d o) {
			// Ascending
			return (Math.abs(x - o.x) > EPS) ? (int) (x - o.x) : (int) (y - o.y);
			// Descending
			// return (Math.abs(x - o.x) > EPS) ? (int) (o.x - x) : (int) (o.y -
			// y);
		}

		public boolean equals(Point_d obj) {
			return (Math.abs(x - obj.x) < EPS) && (Math.abs(y - obj.y) < EPS);
		}
	}

	public static double distance(Point_i p1, Point_i p2) {
		return Math.hypot((double) p1.x - (double) p2.x, (double) p1.y - (double) p2.y);
	}

	public static double distance(Point_d p1, Point_d p2) {
		return Math.hypot(p1.x - p2.x, p1.y - p2.y);
	}

	public static Point_i rotate(Point_i p, double theta) {
		double rad = Math.toRadians(theta); // multiply theta with PI / 180
		return new Point_i((int) (p.x * Math.cos(rad)) - (int) (p.y * Math.sin(rad)),
				(int) (p.x * Math.sin(rad)) - (int) (p.y * Math.cos(rad)));
	}

	public static Point_d rotate(Point_d p, double theta) {
		double rad = Math.toRadians(theta); // multiply theta with PI / 180
		return new Point_d(p.x * Math.cos(rad) - p.y * Math.sin(rad), p.x * Math.sin(rad) - p.y * Math.cos(rad));
	}

	/**
	 * Line
	 */

	public static class LineF {
		// Reprensent Line with Function ax + by + c = 0
		double a;
		double b;
		double c;

		public LineF() {
			a = 0.0;
			b = 0.0;
			c = 0.0;
		}

		public LineF(double a, double b, double c) {
			this.a = a;
			this.b = b;
			this.c = c;
		}
	}

	public static LineF pointsToLineF(Point_i p1, Point_i p2) {
		LineF l = new LineF();

		if (p1.x == p2.x) {
			l.a = 1.0;
			l.b = 0.0;
			l.c = -(double) p1.x;
		} else {
			l.a = -(double) (p1.y - p2.y) / (p1.x - p2.x);
			l.b = 1.0; // IMPORTANT
			l.c = -(double) (l.a * p1.x) - p1.y;
		}

		return l;
	}

	public static LineF pointsToLineF(Point_d p1, Point_d p2) {
		LineF l = new LineF();

		if (Math.abs(p1.x - p2.x) < EPS) {
			l.a = 1.0;
			l.b = 0.0;
			l.c = -p1.x;
		} else {
			l.a = -(p1.y - p2.y) / (p1.x - p2.x);
			l.b = 1.0; // IMPORTANT
			l.c = -(l.a * p1.x) - p1.y;
		}

		return l;
	}

	public static boolean areParallel(LineF l1, LineF l2) {
		return (Math.abs(l1.a - l2.a) < EPS) && (Math.abs(l1.b - l2.b) < EPS);
	}

	public static boolean areSame(LineF l1, LineF l2) {
		return areParallel(l1, l2) && (Math.abs(l1.c - l2.c) < EPS);
	}

	public static Point_d intersect(LineF l1, LineF l2) {
		// Check parallel before use this
		if (areParallel(l1, l2)) {
			return null;
		}

		Point_d p = new Point_d();

		p.x = (l2.b * l1.c - l1.b * l2.c) / (l2.a * l1.b - l1.a * l2.b);
		// Special Case: vertical line to avoid division by 0
		if (Math.abs(l1.b) > EPS) {
			p.y = -(l1.a * p.x + l1.c);
		} else {
			p.y = -(l2.a * p.x + l2.c);
		}

		return p;
	}

	/**
	 * Vec
	 */

	public static class Vec {
		double x;
		double y;

		Vec() {
			x = 0.0;
			y = 0.0;
		}

		Vec(double x, double y) {
			this.x = x;
			this.y = y;
		}
	}

	public static Vec toVec(Point_d a, Point_d b) {
		// Convert 2 points to vector a->b
		return new Vec(b.x - a.x, b.y - a.y);
	}

	public static Vec scale(Vec v, double s) {
		// nonnegative s = [<1 .. 1 .. >1]
		return new Vec(v.x * s, v.y * s);
	}

	public static Point_d translate(Point_d p, Vec v) {
		// translate p according to v
		return new Point_d(p.x + v.x, p.y + v.y);
	}

	public static double dot(Vec a, Vec b) {
		// dot muliply of two vector
		return (a.x * b.x + a.y * b.y);
	}

	public static double norm_sq(Vec v) {
		return (v.x * v.x + v.y * v.y);
	}

	public static double distanceToLine(Point_d p, Point_d a, Point_d b) {
		// Calculate the distance from p to the line defined by a and b
		// formula: c = a + u * ab
		Vec ap = toVec(a, p);
		Vec ab = toVec(a, b);
		double u = dot(ap, ab) / norm_sq(ab);
		Point_d c = translate(p, scale(ab, u));
		return distance(p, c);
	}

	public static Point_d PointToLine(Point_d p, Point_d a, Point_d b) {
		// Get the point p is closest to the Line ab
		// formula: c = a + u * ab
		Vec ap = toVec(a, p);
		Vec ab = toVec(a, b);
		double u = dot(ap, ab) / norm_sq(ab);
		Point_d c = translate(p, scale(ab, u));
		return c;
	}

	public static double distanceToLineSegment(Point_d p, Point_d a, Point_d b) {
		// return the distance from p to the line segment ab defined by point a
		// and b.
		Vec ap = toVec(a, p);
		Vec ab = toVec(a, b);
		double u = dot(ap, ab) / norm_sq(ab);
		if (u < 0.0) {
			return distance(p, a);
		}

		if (u > 1.0) {
			return distance(p, b);
		}

		return distanceToLine(p, a, b);
	}

	public static double angle(Point_d a, Point_d o, Point_d b) {
		// return angle of aob in rad
		// oa * ob = |oa| * |ob| * cos(theta)
		Vec oa = toVec(o, a);
		Vec ob = toVec(o, b);
		return Math.acos(dot(oa, ob) / Math.sqrt(norm_sq(oa) * norm_sq(ob)));
	}

	public static double cross(Vec a, Vec b) {
		return a.x * b.y - a.y * b.x;
	}

	public static boolean ccw(Point_d p, Point_d q, Point_d r) {
		// check if r is on the left side of pq, return true if at left side.
		// If you want to make isConvex work with collinear points, change > to >=
		return cross(toVec(p, q), toVec(p, r)) > 0;
	}

	public static boolean collinear(Point_d p, Point_d q, Point_d r) {
		return Math.abs(cross(toVec(p, q), toVec(p, r))) < EPS;
	}

	/**
	 * Circles
	 */

	public static int insideCircle(Point_i p, Point_i c, int r) {
		// check if point p is in the circle center c and radius r
		int dx = p.x - c.x;
		int dy = p.y - c.y;
		int Euc = dx * dx + dy * dy;
		int rSq = r * r;
		return Euc < rSq ? 0 : Euc == rSq ? 1 : 2; // inside border outside
	}

	public static int insideCircle(Point_d p, Point_d c, double r) {
		// check if point p is in the circle center c and radius r
		double dx = p.x - c.x;
		double dy = p.y - c.y;
		double Euc = dx * dx + dy * dy;
		double rSq = r * r;
		return Euc < rSq ? 0 : Math.abs(Euc - rSq) < EPS ? 1 : 2; // inside
																	// border
																	// outside
	}

	public static double chordLength(double r, double theta) {
		return Math.sqrt(2 * r * r * (1 - Math.cos(Math.toRadians(theta))));
	}

	public static double arcLength(double r, double theta) {
		return (theta / 360) * 2 * PI * r;
	}

	public static double areaSector(double r, double theta) {
		return (theta / 360) * PI * r * r;
	}

	public static Point_d circle2PtsRad(Point_d p1, Point_d p2, double r) {
		// Given two point and the radius
		// Return one of possible circles with p1 p2 on the circle and radius is
		// r
		// Get the other one, just reverse p1 and p2.
		Point_d c = new Point_d();

		double d2 = (p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y);

		double det = r * r / d2 - 0.25;
		if (det < 0.0) {
			return null;
		}
		double h = Math.sqrt(det);
		c.x = (p1.x + p2.x) * 0.5 + (p1.y - p2.y) * h;
		c.y = (p1.y + p2.y) * 0.5 + (p2.x - p1.x) * h;

		return c;
	}

	/**
	 * Triangles
	 */

	public static double perimeter(double a, double b, double c) {
		return a + b + c;
	}

	public static double areaHeronFormula(double a, double b, double c) {
		// semi-perimeter = 0.5 * (a + b + c)
		// area of Triangle with sides a b c. A = Sqrt(s * (s - a) * (s - b) *
		// (s - c))
		double semi = 0.5 * perimeter(a, b, c);
		return Math.sqrt(semi * (semi - a) * (semi - b) * (semi - c));
	}

	public static double rInCircle(double ab, double bc, double ca) {
		// A triangle with area A and semi-perimeter has an inscribed circle
		// (incircle) with radius r = A / s
		return areaHeronFormula(ab, bc, ca) / (0.5 * perimeter(ab, bc, ca));
	}

	public static double rInCircle(Point_d a, Point_d b, Point_d c) {
		return rInCircle(distance(a, b), distance(a, c), distance(b, c));
	}

	public static Point_d inCircle(Point_d p1, Point_d p2, Point_d p3) {
		// return in circle center if exists
		double r = rInCircle(p1, p2, p3);
		if (Math.abs(r) < EPS) {
			return null; // no inCircle center
		}

		LineF l1 = new LineF();
		LineF l2 = new LineF();

		double ratio = distance(p1, p2) / distance(p1, p3);
		Point_d p = translate(p2, scale(toVec(p2, p3), ratio / (1 + ratio)));
		l1 = pointsToLineF(p1, p);

		ratio = distance(p2, p1) / distance(p2, p3);
		p = translate(p1, scale(toVec(p1, p3), ratio / (1 + ratio)));
		l2 = pointsToLineF(p2, p);

		Point_d c = intersect(l1, l2);
		return c;
	}

	public static double rCirCircle(double ab, double bc, double ca) {
		// Return Circumscribed circle radius
		return (ab * bc * ca) / (4 * areaHeronFormula(ab, bc, ca));
	}

	public static double rCirCircle(Point_d a, Point_d b, Point_d c) {
		return rCirCircle(distance(a, b), distance(b, c), distance(c, a));
	}

	/**
	 * Polygon
	 */

	/*
	 * List<Point_d> P = new ArrayList<Point_d>(); P.add(new Point_d(1.0, 2.0));
	 * ... P.add(P[0]) // Important!!!!!!!! loop back
	 */

	public static double perimeter(List<Point_d> P) {
		// Return the perimeter of a Polygon
		double result = 0.0;
		for (int i = 0; i < P.size() - 1; i++) {
			// Remember P[0] = P[n - 1]
			result += distance(P.get(i), P.get(i + 1));
		}
		return result;
	}

	public static double areaPolygon(List<Point_d> P) {
		// Return area of a Polygon
		double result = 0.0;
		double x1 = 0.0;
		double y1 = 0.0;
		double x2 = 0.0;
		double y2 = 0.0;

		for (int i = 0; i < P.size() - 1; i++) {
			x1 = P.get(i).x;
			x2 = P.get(i + 1).x;
			y1 = P.get(i).y;
			y2 = P.get(i + 1).y;
			result += (x1 * y2 - x2 * y1);
		}

		return Math.abs(result) / 2.0;
	}

	public static boolean isConvex(List<Point_d> P) {
		// return true is polygon is convex.
		// Size <= 3 is not convex
		int sz = P.size();
		if (sz <= 3) {
			return false;
		}
		boolean isLeft = ccw(P.get(0), P.get(1), P.get(2)); // remember one
															// result

		for (int i = 1; i < sz - 1; i++) {
			// Compare with the other
			if (ccw(P.get(i), P.get(i + 1), P.get((i + 2) == sz ? 1 : i + 2)) != isLeft) {
				return false;
			}
		}

		return true;
	}

	public static boolean inPolygon(Point_d p, List<Point_d> P) {
		// check if p is in the Polygon P (no matter convex / concave)
		if (P.size() == 0) {
			return false;
		}
		
		double sum = 0; // assume the first vertex is equal to the last vertex
		for (int i = 0; i < P.size() - 1; i++) {
			if (ccw(p, P.get(i), P.get(i + 1))) {
				sum += angle(P.get(i), p, P.get(i + 1)); // left turn ccw
			} else {
				sum -= angle(P.get(i), p, P.get(i + 1)); // right turn ccw
			}
		}

		return Math.abs(Math.abs(sum) - 2 * Math.PI) < EPS;
	}

	public static Point_d lineIntersectSeg(Point_d p, Point_d q, Point_d A, Point_d B) {
		// return point that Line segment p - q intersect with line A - B
		double a = B.y - A.y;
		double b = A.x - B.x;
		double c = B.x * A.y - A.x * B.y;
		double u = Math.abs(a * p.x + b * p.y + c);
		double v = Math.abs(a * q.x + b * q.y + c);

		return new Point_d((p.x * v + q.x * u) / (u + v), (p.y * v + q.y * u) / (u + v));
	}

	public static List<Point_d> cutPolygon(Point_d a, Point_d b, List<Point_d> Q) {
		// Cut polygon Q along the line formed by point a -> point b
		List<Point_d> P = new ArrayList<Point_d>();

		for (int i = 0; i < Q.size(); i++) {
			double left1 = cross(toVec(a, b), toVec(a, Q.get(i)));
			double left2 = 0.0;
			if (i != Q.size() - 1) {
				left2 = cross(toVec(a, b), toVec(a, Q.get(i + 1)));
			}
			if (left1 > -EPS) { // Point Q[i] is on the left of ab
				P.add(Q.get(i));
			}
			if (left1 * left2 < -EPS) { // edge Q[i] Q[i+1] crosses line ab
				P.add(lineIntersectSeg(Q.get(i), Q.get(i + 1), a, b));
			}
		}

		if (!P.isEmpty() && !(P.get(0) == P.get(P.size() - 1))) {
			P.add(P.get(0));
		}

		return P;
	}

	public static Point_d pivot = new Point_d(0, 0);

	public static List<Point_d> CH(List<Point_d> P) {
		// return the convec Hull in a set of Points P
		// the content of P may be reshuffled
		int i = 0;
		int j = 0;
		int n = P.size();

		if (n <= 3) {
			// edge case, the CH is P itself
			if (!(P.get(0) == P.get(n - 1))) {
				P.add(P.get(0));
			}
			return P;
		}

		// First, find P0 = point with lowest Y and if tie: rightmost X;
		int P0 = 0;
		for (i = 1; i < n; i++) {
			if (P.get(i).y < P.get(P0).y || (P.get(i).y == P.get(P0).y && P.get(i).x > P.get(P0).x)) {
				P0 = i;
			}
		}

		Point_d temp = P.get(0);
		P.set(0, P.get(P0));
		P.set(P0, temp);

		// Second, sort points by angle w.r.t. pivot P0
		pivot = P.get(0); // pivot as global
		// Do not sort P[0]
		Collections.sort(P.subList(1, P.size()), new Comparator<Point_d>() {
			@Override
			public int compare(Point_d a, Point_d b) {
				// angle sorting function
				if (collinear(pivot, a, b)) { // special case
					double res = distance(pivot, a) - distance(pivot, b);
					return res > 0 ? 1 : res == 0 ? 0 : -1; // check which one
															// is close
				}
				double d1x = a.x - pivot.x;
				double d1y = a.y - pivot.y;
				double d2x = b.x - pivot.x;
				double d2y = b.y - pivot.y;

				double r = Math.atan2(d1y, d1x) - Math.atan2(d2y, d2x);
				return r > 0 ? 1 : r == 0 ? 0 : -1; // Compare
													// two
													// angles
			}
		});
		
		// Third  the ccw test
		List<Point_d> S = new ArrayList<Point_d>();
		
		S.add(P.get(n - 1));
		S.add(P.get(0));
		S.add(P.get(1));
		i = 2;
		while (i < n) { //N must be >= 3 for this method to work
			j = S.size() - 1;
			if(ccw(S.get(j - 1), S.get(j), P.get(i))) {
				S.add(P.get(i++)); // left turn accept
			} else {
				S.remove(S.size() - 1); // pop the top of S until we have a left turn
			}
		}
		
		return S;
	}
	
	public static void main(String[] args) {
		Point_d a = new Point_d(0, 0);
		Point_d b = new Point_d(1, 0);
		Point_d c = new Point_d(2.1, 0.9);
		Point_d d = new Point_d(3, 2);
		Point_d e = new Point_d(3, 3);
		Point_d f = new Point_d(1, 3);
		
		List<Point_d> P = new ArrayList<Point_d>();
		P.add(a);
		P.add(b);
		P.add(c);
		P.add(d);
		P.add(e);
		P.add(f);
		P.add(a);
		
		for (int i = 0; i < P.size(); i++) {
			System.out.println(P.get(i).x + "-" + P.get(i).y);
		}
		
		System.out.println("Perimeter: " + perimeter(P));
		System.out.println("Area: " + areaPolygon(P));
		System.out.println("isConvex: " + isConvex(P));
		
		Point_d g = new Point_d(1, 1);
		Point_d h = new Point_d(2, 2);
		
		List<Point_d> Q = new ArrayList<Point_d>();
		Q.add(a);
		Q.add(b);
		Q.add(c);
		Q.add(d);
		Q.add(e);
		Q.add(f);
		Q.add(g);
		Q.add(h);
		
		List<Point_d> res = CH(Q);
		
		for (int i = 0; i < res.size(); i++) {
			System.out.println(res.get(i).x + "-" + res.get(i).y);
		}
	}
}
