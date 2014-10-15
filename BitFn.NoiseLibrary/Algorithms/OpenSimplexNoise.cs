namespace BitFn.NoiseLibrary.Algorithms
{
	/// <summary>
	///     A port of Kurt Spencer's OpenSimplexNoise.java. Computes noise using surflet summation on a Simplectic Honeycomb.
	/// </summary>
	/// <remarks>
	///     Original Java implementation released under public domain:
	///     https://gist.github.com/KdotJPG/b1270127455a94ac5d19/b6fb5128cf983dffa74e727dbcdebaa23dfce5f4
	/// </remarks>
	public class OpenSimplexNoise : INoiseService2D, INoiseService3D, INoiseService4D
	{
		/// <summary>
		///     The squish constant to use for two dimensional noise.
		/// </summary>
		/// <remarks>
		///     <code>(sqrt(2+1)-1)/2 = 0.366025403784438646763723170752936183471402626905190314027903</code>
		/// </remarks>
		private const double Squish2D = 0.366025403784438646763723170752936183471402626905190314027903;

		/// <summary>
		///     The squish constant to use for three dimensional noise.
		/// </summary>
		/// <remarks>
		///     <code>(sqrt(3+1)-1)/3 = 1/3</code>
		/// </remarks>
		private const double Squish3D = 1D/3;

		/// <summary>
		///     The squish constant to use for four dimensional noise.
		/// </summary>
		/// <remarks>
		///     <code>(sqrt(4+1)-1)/4 = 0.309016994374947424102293417182819058860154589902881431067724</code>
		/// </remarks>
		private const double Squish4D = 0.309016994374947424102293417182819058860154589902881431067724;

		/// <summary>
		///     The strech constant to use for two dimensional noise.
		/// </summary>
		/// <remarks>
		///     <code>(1/sqrt(2+1)-1)/2 = -0.21132486540518711774542560974902127217619912436493656199069</code>
		/// </remarks>
		private const double Stretch2D = -0.21132486540518711774542560974902127217619912436493656199069;

		/// <summary>
		///     The strech constant to use for three dimensional noise.
		/// </summary>
		/// <remarks>
		///     <code>(1/sqrt(3+1)-1)/3 = -1/6</code>
		/// </remarks>
		private const double Stretch3D = -1D/6;

		/// <summary>
		///     The strech constant to use for four dimensional noise.
		/// </summary>
		/// <remarks>
		///     <code>(1/sqrt(4+1)-1)/4 = -0.13819660112501051517954131656343618822796908201942371378645</code>
		/// </remarks>
		private const double Stretch4D = -0.13819660112501051517954131656343618822796908201942371378645;

		/// <summary>
		///     The normalization scalar to use for two dimensional noise.
		/// </summary>
		private const double Normalization2D = 47;

		/// <summary>
		///     The normalization scalar to use for three dimensional noise.
		/// </summary>
		private const double Normalization3D = 103;

		/// <summary>
		///     The normalization scalar to use for four dimensional noise.
		/// </summary>
		private const double Normalization4D = 30;

		/// <summary>
		///     The default seed to use for the parameterless constructor.
		/// </summary>
		private const long DefaultSeed = 0;

		private const int PermutationLength = 256;
		private const long PermutationConstant1 = 6364136223846793005L;
		private const long PermutationConstant2 = 1442695040888963407L;

		/// <summary>
		///     Gradients for 2D. They approximate the directions to the vertices of an octagon from the center.
		/// </summary>
		private static readonly sbyte[] Gradients2D =
		{
			5, 2, 2, 5, -5, 2, -2, 5, 5, -2, 2, -5, -5, -2, -2, -5
		};

		/// <summary>
		///     Gradients for 3D. They approximate the directions to the vertices of a rhombicuboctahedron from the center, skewed
		///     so that the triangular and square facets can be inscribed inside circles of the same radius.
		/// </summary>
		private static readonly sbyte[] Gradients3D =
		{
			-11, 4, 4, -4, 11, 4, -4, 4, 11, 11, 4, 4, 4, 11, 4, 4, 4, 11, -11, -4, 4, -4, -11, 4, -4, -4, 11, 11, -4, 4, 4, -11,
			4, 4, -4, 11, -11, 4, -4, -4, 11, -4, -4, 4, -11, 11, 4, -4, 4, 11, -4, 4, 4, -11, -11, -4, -4, -4, -11, -4, -4, -4,
			-11, 11, -4, -4, 4, -11, -4, 4, -4, -11
		};


		/// <summary>
		///     Gradients for 4D. They approximate the directions to the vertices of a disprismatotesseractihexadecachoron from the
		///     center, skewed so that the tetrahedral and cubic facets can be inscribed inside spheres of the same radius.
		/// </summary>
		private static readonly sbyte[] Gradients4D =
		{
			3, 1, 1, 1, 1, 3, 1, 1, 1, 1, 3, 1, 1, 1, 1, 3, -3, 1, 1, 1, -1, 3, 1, 1, -1, 1, 3, 1, -1, 1, 1, 3, 3, -1, 1, 1, 1,
			-3, 1, 1, 1, -1, 3, 1, 1, -1, 1, 3, -3, -1, 1, 1, -1, -3, 1, 1, -1, -1, 3, 1, -1, -1, 1, 3, 3, 1, -1, 1, 1, 3, -1, 1,
			1, 1, -3, 1, 1, 1, -1, 3, -3, 1, -1, 1, -1, 3, -1, 1, -1, 1, -3, 1, -1, 1, -1, 3, 3, -1, -1, 1, 1, -3, -1, 1, 1, -1,
			-3, 1, 1, -1, -1, 3, -3, -1, -1, 1, -1, -3, -1, 1, -1, -1, -3, 1, -1, -1, -1, 3, 3, 1, 1, -1, 1, 3, 1, -1, 1, 1, 3,
			-1, 1, 1, 1, -3, -3, 1, 1, -1, -1, 3, 1, -1, -1, 1, 3, -1, -1, 1, 1, -3, 3, -1, 1, -1, 1, -3, 1, -1, 1, -1, 3, -1, 1,
			-1, 1, -3, -3, -1, 1, -1, -1, -3, 1, -1, -1, -1, 3, -1, -1, -1, 1, -3, 3, 1, -1, -1, 1, 3, -1, -1, 1, 1, -3, -1, 1, 1,
			-1, -3, -3, 1, -1, -1, -1, 3, -1, -1, -1, 1, -3, -1, -1, 1, -1, -3, 3, -1, -1, -1, 1, -3, -1, -1, 1, -1, -3, -1, 1,
			-1, -1, -3, -3, -1, -1, -1, -1, -3, -1, -1, -1, -1, -3, -1, -1, -1, -1, -3
		};

		private readonly byte[] _permutation;
		private readonly byte[] _permutationGradIndex3D;

		public OpenSimplexNoise()
			: this(DefaultSeed)
		{
		}

		public OpenSimplexNoise(byte[] perm)
		{
			_permutation = perm;
			_permutationGradIndex3D = new byte[PermutationLength];

			for (int i = 0; i < PermutationLength; i++)
			{
				// Since 3D has 24 gradients, simple bitmask won't work, so precompute modulo array.
				_permutationGradIndex3D[i] = (byte) ((perm[i]%(Gradients3D.Length/3))*3);
			}
		}

		/// <summary>
		///     Initializes the class using a permutation array generated from a 64-bit seed. Generates a proper permutation (i.e.
		///     doesn't merely perform N successive pair swaps on a base array) Uses a simple 64-bit LCG.
		/// </summary>
		/// <param name="seed"></param>
		public OpenSimplexNoise(long seed)
		{
			_permutation = new byte[PermutationLength];
			_permutationGradIndex3D = new byte[PermutationLength];
			var source = new byte[PermutationLength];
			for (int i = 0; i < PermutationLength; i++)
				source[i] = (byte) i;
			seed = seed*PermutationConstant1 + PermutationConstant2;
			seed = seed*PermutationConstant1 + PermutationConstant2;
			seed = seed*PermutationConstant1 + PermutationConstant2;
			for (int i = PermutationLength - 1; i >= 0; i--)
			{
				seed = seed*PermutationConstant1 + PermutationConstant2;
				var r = (int) ((seed + 31)%(i + 1));
				if (r < 0)
					r += (i + 1);
				_permutation[i] = source[r];
				_permutationGradIndex3D[i] = (byte) ((_permutation[i]%(Gradients3D.Length/3))*3);
				source[r] = source[i];
			}
		}

		public double ComputeNoise(double x, double y)
		{
			// Place input coordinates onto grid.
			double stretchOffset = (x + y)*Stretch2D;
			double xs = x + stretchOffset;
			double ys = y + stretchOffset;

			// Floor to get grid coordinates of rhombus (stretched square) super-cell origin.
			int xsb = FastFloor(xs);
			int ysb = FastFloor(ys);

			// Skew out to get actual coordinates of rhombus origin. We'll need these later.
			double squishOffset = (xsb + ysb)*Squish2D;
			double xb = xsb + squishOffset;
			double yb = ysb + squishOffset;

			// Compute grid coordinates relative to rhombus origin.
			double xins = xs - xsb;
			double yins = ys - ysb;

			// Sum those together to get a value that determines which region we're in.
			double inSum = xins + yins;

			// Positions relative to origin point.
			double dx0 = x - xb;
			double dy0 = y - yb;

			// We'll be defining these inside the next block and using them afterwards.
			double dxExt, dyExt;
			int xsvExt, ysvExt;

			double value = 0;

			// Contribution (1,0)
			double dx1 = dx0 - 1 - Squish2D;
			double dy1 = dy0 - 0 - Squish2D;
			double attn1 = 2 - dx1*dx1 - dy1*dy1;
			if (attn1 > 0)
			{
				attn1 *= attn1;
				value += attn1*attn1*Extrapolate(xsb + 1, ysb + 0, dx1, dy1);
			}

			// Contribution (0,1)
			double dx2 = dx0 - 0 - Squish2D;
			double dy2 = dy0 - 1 - Squish2D;
			double attn2 = 2 - dx2*dx2 - dy2*dy2;
			if (attn2 > 0)
			{
				attn2 *= attn2;
				value += attn2*attn2*Extrapolate(xsb + 0, ysb + 1, dx2, dy2);
			}

			if (inSum <= 1)
			{
				// We're inside the triangle (2-Simplex) at (0,0)
				double zins = 1 - inSum;
				if (zins > xins || zins > yins)
				{
					// (0,0) is one of the closest two triangular vertices
					if (xins > yins)
					{
						xsvExt = xsb + 1;
						ysvExt = ysb - 1;
						dxExt = dx0 - 1;
						dyExt = dy0 + 1;
					}
					else
					{
						xsvExt = xsb - 1;
						ysvExt = ysb + 1;
						dxExt = dx0 + 1;
						dyExt = dy0 - 1;
					}
				}
				else
				{
					// (1,0) and (0,1) are the closest two vertices.
					xsvExt = xsb + 1;
					ysvExt = ysb + 1;
					dxExt = dx0 - 1 - 2*Squish2D;
					dyExt = dy0 - 1 - 2*Squish2D;
				}
			}
			else
			{
				// We're inside the triangle (2-Simplex) at (1,1)
				double zins = 2 - inSum;
				if (zins < xins || zins < yins)
				{
					// (0,0) is one of the closest two triangular vertices
					if (xins > yins)
					{
						xsvExt = xsb + 2;
						ysvExt = ysb + 0;
						dxExt = dx0 - 2 - 2*Squish2D;
						dyExt = dy0 + 0 - 2*Squish2D;
					}
					else
					{
						xsvExt = xsb + 0;
						ysvExt = ysb + 2;
						dxExt = dx0 + 0 - 2*Squish2D;
						dyExt = dy0 - 2 - 2*Squish2D;
					}
				}
				else
				{
					// (1,0) and (0,1) are the closest two vertices.
					dxExt = dx0;
					dyExt = dy0;
					xsvExt = xsb;
					ysvExt = ysb;
				}
				xsb += 1;
				ysb += 1;
				dx0 = dx0 - 1 - 2*Squish2D;
				dy0 = dy0 - 1 - 2*Squish2D;
			}

			// Contribution (0,0) or (1,1)
			double attn0 = 2 - dx0*dx0 - dy0*dy0;
			if (attn0 > 0)
			{
				attn0 *= attn0;
				value += attn0*attn0*Extrapolate(xsb, ysb, dx0, dy0);
			}

			// Extra Vertex
			double attnExt = 2 - dxExt*dxExt - dyExt*dyExt;
			if (attnExt > 0)
			{
				attnExt *= attnExt;
				value += attnExt*attnExt*Extrapolate(xsvExt, ysvExt, dxExt, dyExt);
			}

			return value/Normalization2D;
		}

		// 3D OpenSimplex (Simplectic) Noise.
		public double ComputeNoise(double x, double y, double z)
		{
			// Place input coordinates on simplectic honeycomb.
			double stretchOffset = (x + y + z)*Stretch3D;
			double xs = x + stretchOffset;
			double ys = y + stretchOffset;
			double zs = z + stretchOffset;

			// Floor to get simplectic honeycomb coordinates of rhombohedron (stretched cube) super-cell origin.
			int xsb = FastFloor(xs);
			int ysb = FastFloor(ys);
			int zsb = FastFloor(zs);

			// Skew out to get actual coordinates of rhombohedron origin. We'll need these later.
			double squishOffset = (xsb + ysb + zsb)*Squish3D;
			double xb = xsb + squishOffset;
			double yb = ysb + squishOffset;
			double zb = zsb + squishOffset;

			// Compute simplectic honeycomb coordinates relative to rhombohedral origin.
			double xins = xs - xsb;
			double yins = ys - ysb;
			double zins = zs - zsb;

			// Sum those together to get a value that determines which region we're in.
			double inSum = xins + yins + zins;

			// Positions relative to origin point.
			double dx0 = x - xb;
			double dy0 = y - yb;
			double dz0 = z - zb;

			// We'll be defining these inside the next block and using them afterwards.
			double dxExt0, dyExt0, dzExt0;
			double dxExt1, dyExt1, dzExt1;
			int xsvExt0, ysvExt0, zsvExt0;
			int xsvExt1, ysvExt1, zsvExt1;

			double value = 0;
			if (inSum <= 1)
			{
				// We're inside the tetrahedron (3-Simplex) at (0,0,0)

				// Determine which two of (0,0,1), (0,1,0), (1,0,0) are closest.
				sbyte aPoint = 0x01;
				double aScore = xins;
				sbyte bPoint = 0x02;
				double bScore = yins;
				if (aScore >= bScore && zins > bScore)
				{
					bScore = zins;
					bPoint = 0x04;
				}
				else if (aScore < bScore && zins > aScore)
				{
					aScore = zins;
					aPoint = 0x04;
				}

				// Now we determine the two lattice points not part of the tetrahedron that may contribute.
				// This depends on the closest two tetrahedral vertices, including (0,0,0)
				double wins = 1 - inSum;
				if (wins > aScore || wins > bScore)
				{
					// (0,0,0) is one of the closest two tetrahedral vertices.
					sbyte c = (bScore > aScore ? bPoint : aPoint); // Our other closest vertex is the closest out of a and b.

					if ((c & 0x01) == 0)
					{
						xsvExt0 = xsb - 1;
						xsvExt1 = xsb;
						dxExt0 = dx0 + 1;
						dxExt1 = dx0;
					}
					else
					{
						xsvExt0 = xsvExt1 = xsb + 1;
						dxExt0 = dxExt1 = dx0 - 1;
					}

					if ((c & 0x02) == 0)
					{
						ysvExt0 = ysvExt1 = ysb;
						dyExt0 = dyExt1 = dy0;
						if ((c & 0x01) == 0)
						{
							ysvExt1 -= 1;
							dyExt1 += 1;
						}
						else
						{
							ysvExt0 -= 1;
							dyExt0 += 1;
						}
					}
					else
					{
						ysvExt0 = ysvExt1 = ysb + 1;
						dyExt0 = dyExt1 = dy0 - 1;
					}

					if ((c & 0x04) == 0)
					{
						zsvExt0 = zsb;
						zsvExt1 = zsb - 1;
						dzExt0 = dz0;
						dzExt1 = dz0 + 1;
					}
					else
					{
						zsvExt0 = zsvExt1 = zsb + 1;
						dzExt0 = dzExt1 = dz0 - 1;
					}
				}
				else
				{
					// (0,0,0) is not one of the closest two tetrahedral vertices.
					var c = (sbyte) (aPoint | bPoint); // Our two extra vertices are determined by the closest two.

					if ((c & 0x01) == 0)
					{
						xsvExt0 = xsb;
						xsvExt1 = xsb - 1;
						dxExt0 = dx0 - 2*Squish3D;
						dxExt1 = dx0 + 1 - Squish3D;
					}
					else
					{
						xsvExt0 = xsvExt1 = xsb + 1;
						dxExt0 = dx0 - 1 - 2*Squish3D;
						dxExt1 = dx0 - 1 - Squish3D;
					}

					if ((c & 0x02) == 0)
					{
						ysvExt0 = ysb;
						ysvExt1 = ysb - 1;
						dyExt0 = dy0 - 2*Squish3D;
						dyExt1 = dy0 + 1 - Squish3D;
					}
					else
					{
						ysvExt0 = ysvExt1 = ysb + 1;
						dyExt0 = dy0 - 1 - 2*Squish3D;
						dyExt1 = dy0 - 1 - Squish3D;
					}

					if ((c & 0x04) == 0)
					{
						zsvExt0 = zsb;
						zsvExt1 = zsb - 1;
						dzExt0 = dz0 - 2*Squish3D;
						dzExt1 = dz0 + 1 - Squish3D;
					}
					else
					{
						zsvExt0 = zsvExt1 = zsb + 1;
						dzExt0 = dz0 - 1 - 2*Squish3D;
						dzExt1 = dz0 - 1 - Squish3D;
					}
				}

				// Contribution (0,0,0)
				double attn0 = 2 - dx0*dx0 - dy0*dy0 - dz0*dz0;
				if (attn0 > 0)
				{
					attn0 *= attn0;
					value += attn0*attn0*Extrapolate(xsb + 0, ysb + 0, zsb + 0, dx0, dy0, dz0);
				}

				// Contribution (1,0,0)
				double dx1 = dx0 - 1 - Squish3D;
				double dy1 = dy0 - 0 - Squish3D;
				double dz1 = dz0 - 0 - Squish3D;
				double attn1 = 2 - dx1*dx1 - dy1*dy1 - dz1*dz1;
				if (attn1 > 0)
				{
					attn1 *= attn1;
					value += attn1*attn1*Extrapolate(xsb + 1, ysb + 0, zsb + 0, dx1, dy1, dz1);
				}

				// Contribution (0,1,0)
				double dx2 = dx0 - 0 - Squish3D;
				double dy2 = dy0 - 1 - Squish3D;
				double dz2 = dz1;
				double attn2 = 2 - dx2*dx2 - dy2*dy2 - dz2*dz2;
				if (attn2 > 0)
				{
					attn2 *= attn2;
					value += attn2*attn2*Extrapolate(xsb + 0, ysb + 1, zsb + 0, dx2, dy2, dz2);
				}

				// Contribution (0,0,1)
				double dx3 = dx2;
				double dy3 = dy1;
				double dz3 = dz0 - 1 - Squish3D;
				double attn3 = 2 - dx3*dx3 - dy3*dy3 - dz3*dz3;
				if (attn3 > 0)
				{
					attn3 *= attn3;
					value += attn3*attn3*Extrapolate(xsb + 0, ysb + 0, zsb + 1, dx3, dy3, dz3);
				}
			}
			else if (inSum >= 2)
			{
				// We're inside the tetrahedron (3-Simplex) at (1,1,1)

				// Determine which two tetrahedral vertices are the closest, out of (1,1,0), (1,0,1), (0,1,1) but not (1,1,1).
				sbyte aPoint = 0x06;
				double aScore = xins;
				sbyte bPoint = 0x05;
				double bScore = yins;
				if (aScore <= bScore && zins < bScore)
				{
					bScore = zins;
					bPoint = 0x03;
				}
				else if (aScore > bScore && zins < aScore)
				{
					aScore = zins;
					aPoint = 0x03;
				}

				// Now we determine the two lattice points not part of the tetrahedron that may contribute.
				// This depends on the closest two tetrahedral vertices, including (1,1,1)
				double wins = 3 - inSum;
				if (wins < aScore || wins < bScore)
				{
					// (1,1,1) is one of the closest two tetrahedral vertices.
					sbyte c = (bScore < aScore ? bPoint : aPoint); // Our other closest vertex is the closest out of a and b.

					if ((c & 0x01) != 0)
					{
						xsvExt0 = xsb + 2;
						xsvExt1 = xsb + 1;
						dxExt0 = dx0 - 2 - 3*Squish3D;
						dxExt1 = dx0 - 1 - 3*Squish3D;
					}
					else
					{
						xsvExt0 = xsvExt1 = xsb;
						dxExt0 = dxExt1 = dx0 - 3*Squish3D;
					}

					if ((c & 0x02) != 0)
					{
						ysvExt0 = ysvExt1 = ysb + 1;
						dyExt0 = dyExt1 = dy0 - 1 - 3*Squish3D;
						if ((c & 0x01) != 0)
						{
							ysvExt1 += 1;
							dyExt1 -= 1;
						}
						else
						{
							ysvExt0 += 1;
							dyExt0 -= 1;
						}
					}
					else
					{
						ysvExt0 = ysvExt1 = ysb;
						dyExt0 = dyExt1 = dy0 - 3*Squish3D;
					}

					if ((c & 0x04) != 0)
					{
						zsvExt0 = zsb + 1;
						zsvExt1 = zsb + 2;
						dzExt0 = dz0 - 1 - 3*Squish3D;
						dzExt1 = dz0 - 2 - 3*Squish3D;
					}
					else
					{
						zsvExt0 = zsvExt1 = zsb;
						dzExt0 = dzExt1 = dz0 - 3*Squish3D;
					}
				}
				else
				{
					// (1,1,1) is not one of the closest two tetrahedral vertices.
					var c = (sbyte) (aPoint & bPoint); // Our two extra vertices are determined by the closest two.

					if ((c & 0x01) != 0)
					{
						xsvExt0 = xsb + 1;
						xsvExt1 = xsb + 2;
						dxExt0 = dx0 - 1 - Squish3D;
						dxExt1 = dx0 - 2 - 2*Squish3D;
					}
					else
					{
						xsvExt0 = xsvExt1 = xsb;
						dxExt0 = dx0 - Squish3D;
						dxExt1 = dx0 - 2*Squish3D;
					}

					if ((c & 0x02) != 0)
					{
						ysvExt0 = ysb + 1;
						ysvExt1 = ysb + 2;
						dyExt0 = dy0 - 1 - Squish3D;
						dyExt1 = dy0 - 2 - 2*Squish3D;
					}
					else
					{
						ysvExt0 = ysvExt1 = ysb;
						dyExt0 = dy0 - Squish3D;
						dyExt1 = dy0 - 2*Squish3D;
					}

					if ((c & 0x04) != 0)
					{
						zsvExt0 = zsb + 1;
						zsvExt1 = zsb + 2;
						dzExt0 = dz0 - 1 - Squish3D;
						dzExt1 = dz0 - 2 - 2*Squish3D;
					}
					else
					{
						zsvExt0 = zsvExt1 = zsb;
						dzExt0 = dz0 - Squish3D;
						dzExt1 = dz0 - 2*Squish3D;
					}
				}

				// Contribution (1,1,0)
				double dx3 = dx0 - 1 - 2*Squish3D;
				double dy3 = dy0 - 1 - 2*Squish3D;
				double dz3 = dz0 - 0 - 2*Squish3D;
				double attn3 = 2 - dx3*dx3 - dy3*dy3 - dz3*dz3;
				if (attn3 > 0)
				{
					attn3 *= attn3;
					value += attn3*attn3*Extrapolate(xsb + 1, ysb + 1, zsb + 0, dx3, dy3, dz3);
				}

				// Contribution (1,0,1)
				double dx2 = dx3;
				double dy2 = dy0 - 0 - 2*Squish3D;
				double dz2 = dz0 - 1 - 2*Squish3D;
				double attn2 = 2 - dx2*dx2 - dy2*dy2 - dz2*dz2;
				if (attn2 > 0)
				{
					attn2 *= attn2;
					value += attn2*attn2*Extrapolate(xsb + 1, ysb + 0, zsb + 1, dx2, dy2, dz2);
				}

				// Contribution (0,1,1)
				double dx1 = dx0 - 0 - 2*Squish3D;
				double dy1 = dy3;
				double dz1 = dz2;
				double attn1 = 2 - dx1*dx1 - dy1*dy1 - dz1*dz1;
				if (attn1 > 0)
				{
					attn1 *= attn1;
					value += attn1*attn1*Extrapolate(xsb + 0, ysb + 1, zsb + 1, dx1, dy1, dz1);
				}

				// Contribution (1,1,1)
				dx0 = dx0 - 1 - 3*Squish3D;
				dy0 = dy0 - 1 - 3*Squish3D;
				dz0 = dz0 - 1 - 3*Squish3D;
				double attn0 = 2 - dx0*dx0 - dy0*dy0 - dz0*dz0;
				if (attn0 > 0)
				{
					attn0 *= attn0;
					value += attn0*attn0*Extrapolate(xsb + 1, ysb + 1, zsb + 1, dx0, dy0, dz0);
				}
			}
			else
			{
				// We're inside the octahedron (Rectified 3-Simplex) in between.
				double aScore;
				sbyte aPoint;
				bool aIsFurtherSide;
				double bScore;
				sbyte bPoint;
				bool bIsFurtherSide;

				// Decide between point (0,0,1) and (1,1,0) as closest
				double p1 = xins + yins;
				if (p1 > 1)
				{
					aScore = p1 - 1;
					aPoint = 0x03;
					aIsFurtherSide = true;
				}
				else
				{
					aScore = 1 - p1;
					aPoint = 0x04;
					aIsFurtherSide = false;
				}

				// Decide between point (0,1,0) and (1,0,1) as closest
				double p2 = xins + zins;
				if (p2 > 1)
				{
					bScore = p2 - 1;
					bPoint = 0x05;
					bIsFurtherSide = true;
				}
				else
				{
					bScore = 1 - p2;
					bPoint = 0x02;
					bIsFurtherSide = false;
				}

				// The closest out of the two (1,0,0) and (0,1,1) will replace the furthest out of the two decided above, if closer.
				double p3 = yins + zins;
				if (p3 > 1)
				{
					double score = p3 - 1;
					if (aScore <= bScore && aScore < score)
					{
						// aScore = score;
						aPoint = 0x06;
						aIsFurtherSide = true;
					}
					else if (aScore > bScore && bScore < score)
					{
						// bScore = score;
						bPoint = 0x06;
						bIsFurtherSide = true;
					}
				}
				else
				{
					double score = 1 - p3;
					if (aScore <= bScore && aScore < score)
					{
						// aScore = score;
						aPoint = 0x01;
						aIsFurtherSide = false;
					}
					else if (aScore > bScore && bScore < score)
					{
						// bScore = score;
						bPoint = 0x01;
						bIsFurtherSide = false;
					}
				}

				// Where each of the two closest points are determines how the extra two vertices are calculated.
				if (aIsFurtherSide == bIsFurtherSide)
				{
					if (aIsFurtherSide)
					{
						// Both closest points on (1,1,1) side

						// One of the two extra points is (1,1,1)
						dxExt0 = dx0 - 1 - 3*Squish3D;
						dyExt0 = dy0 - 1 - 3*Squish3D;
						dzExt0 = dz0 - 1 - 3*Squish3D;
						xsvExt0 = xsb + 1;
						ysvExt0 = ysb + 1;
						zsvExt0 = zsb + 1;

						// Other extra point is based on the shared axis.
						var c = (sbyte) (aPoint & bPoint);
						if ((c & 0x01) != 0)
						{
							dxExt1 = dx0 - 2 - 2*Squish3D;
							dyExt1 = dy0 - 2*Squish3D;
							dzExt1 = dz0 - 2*Squish3D;
							xsvExt1 = xsb + 2;
							ysvExt1 = ysb;
							zsvExt1 = zsb;
						}
						else if ((c & 0x02) != 0)
						{
							dxExt1 = dx0 - 2*Squish3D;
							dyExt1 = dy0 - 2 - 2*Squish3D;
							dzExt1 = dz0 - 2*Squish3D;
							xsvExt1 = xsb;
							ysvExt1 = ysb + 2;
							zsvExt1 = zsb;
						}
						else
						{
							dxExt1 = dx0 - 2*Squish3D;
							dyExt1 = dy0 - 2*Squish3D;
							dzExt1 = dz0 - 2 - 2*Squish3D;
							xsvExt1 = xsb;
							ysvExt1 = ysb;
							zsvExt1 = zsb + 2;
						}
					}
					else
					{
						// Both closest points on (0,0,0) side

						// One of the two extra points is (0,0,0)
						dxExt0 = dx0;
						dyExt0 = dy0;
						dzExt0 = dz0;
						xsvExt0 = xsb;
						ysvExt0 = ysb;
						zsvExt0 = zsb;

						// Other extra point is based on the omitted axis.
						var c = (sbyte) (aPoint | bPoint);
						if ((c & 0x01) == 0)
						{
							dxExt1 = dx0 + 1 - Squish3D;
							dyExt1 = dy0 - 1 - Squish3D;
							dzExt1 = dz0 - 1 - Squish3D;
							xsvExt1 = xsb - 1;
							ysvExt1 = ysb + 1;
							zsvExt1 = zsb + 1;
						}
						else if ((c & 0x02) == 0)
						{
							dxExt1 = dx0 - 1 - Squish3D;
							dyExt1 = dy0 + 1 - Squish3D;
							dzExt1 = dz0 - 1 - Squish3D;
							xsvExt1 = xsb + 1;
							ysvExt1 = ysb - 1;
							zsvExt1 = zsb + 1;
						}
						else
						{
							dxExt1 = dx0 - 1 - Squish3D;
							dyExt1 = dy0 - 1 - Squish3D;
							dzExt1 = dz0 + 1 - Squish3D;
							xsvExt1 = xsb + 1;
							ysvExt1 = ysb + 1;
							zsvExt1 = zsb - 1;
						}
					}
				}
				else
				{
					// One point on (0,0,0) side, one point on (1,1,1) side
					sbyte c1, c2;
					if (aIsFurtherSide)
					{
						c1 = aPoint;
						c2 = bPoint;
					}
					else
					{
						c1 = bPoint;
						c2 = aPoint;
					}

					// One contribution is a permutation of (1,1,-1)
					if ((c1 & 0x01) == 0)
					{
						dxExt0 = dx0 + 1 - Squish3D;
						dyExt0 = dy0 - 1 - Squish3D;
						dzExt0 = dz0 - 1 - Squish3D;
						xsvExt0 = xsb - 1;
						ysvExt0 = ysb + 1;
						zsvExt0 = zsb + 1;
					}
					else if ((c1 & 0x02) == 0)
					{
						dxExt0 = dx0 - 1 - Squish3D;
						dyExt0 = dy0 + 1 - Squish3D;
						dzExt0 = dz0 - 1 - Squish3D;
						xsvExt0 = xsb + 1;
						ysvExt0 = ysb - 1;
						zsvExt0 = zsb + 1;
					}
					else
					{
						dxExt0 = dx0 - 1 - Squish3D;
						dyExt0 = dy0 - 1 - Squish3D;
						dzExt0 = dz0 + 1 - Squish3D;
						xsvExt0 = xsb + 1;
						ysvExt0 = ysb + 1;
						zsvExt0 = zsb - 1;
					}

					// One contribution is a permutation of (0,0,2)
					dxExt1 = dx0 - 2*Squish3D;
					dyExt1 = dy0 - 2*Squish3D;
					dzExt1 = dz0 - 2*Squish3D;
					xsvExt1 = xsb;
					ysvExt1 = ysb;
					zsvExt1 = zsb;
					if ((c2 & 0x01) != 0)
					{
						dxExt1 -= 2;
						xsvExt1 += 2;
					}
					else if ((c2 & 0x02) != 0)
					{
						dyExt1 -= 2;
						ysvExt1 += 2;
					}
					else
					{
						dzExt1 -= 2;
						zsvExt1 += 2;
					}
				}

				// Contribution (1,0,0)
				double dx1 = dx0 - 1 - Squish3D;
				double dy1 = dy0 - 0 - Squish3D;
				double dz1 = dz0 - 0 - Squish3D;
				double attn1 = 2 - dx1*dx1 - dy1*dy1 - dz1*dz1;
				if (attn1 > 0)
				{
					attn1 *= attn1;
					value += attn1*attn1*Extrapolate(xsb + 1, ysb + 0, zsb + 0, dx1, dy1, dz1);
				}

				// Contribution (0,1,0)
				double dx2 = dx0 - 0 - Squish3D;
				double dy2 = dy0 - 1 - Squish3D;
				double dz2 = dz1;
				double attn2 = 2 - dx2*dx2 - dy2*dy2 - dz2*dz2;
				if (attn2 > 0)
				{
					attn2 *= attn2;
					value += attn2*attn2*Extrapolate(xsb + 0, ysb + 1, zsb + 0, dx2, dy2, dz2);
				}

				// Contribution (0,0,1)
				double dx3 = dx2;
				double dy3 = dy1;
				double dz3 = dz0 - 1 - Squish3D;
				double attn3 = 2 - dx3*dx3 - dy3*dy3 - dz3*dz3;
				if (attn3 > 0)
				{
					attn3 *= attn3;
					value += attn3*attn3*Extrapolate(xsb + 0, ysb + 0, zsb + 1, dx3, dy3, dz3);
				}

				// Contribution (1,1,0)
				double dx4 = dx0 - 1 - 2*Squish3D;
				double dy4 = dy0 - 1 - 2*Squish3D;
				double dz4 = dz0 - 0 - 2*Squish3D;
				double attn4 = 2 - dx4*dx4 - dy4*dy4 - dz4*dz4;
				if (attn4 > 0)
				{
					attn4 *= attn4;
					value += attn4*attn4*Extrapolate(xsb + 1, ysb + 1, zsb + 0, dx4, dy4, dz4);
				}

				// Contribution (1,0,1)
				double dx5 = dx4;
				double dy5 = dy0 - 0 - 2*Squish3D;
				double dz5 = dz0 - 1 - 2*Squish3D;
				double attn5 = 2 - dx5*dx5 - dy5*dy5 - dz5*dz5;
				if (attn5 > 0)
				{
					attn5 *= attn5;
					value += attn5*attn5*Extrapolate(xsb + 1, ysb + 0, zsb + 1, dx5, dy5, dz5);
				}

				// Contribution (0,1,1)
				double dx6 = dx0 - 0 - 2*Squish3D;
				double dy6 = dy4;
				double dz6 = dz5;
				double attn6 = 2 - dx6*dx6 - dy6*dy6 - dz6*dz6;
				if (attn6 > 0)
				{
					attn6 *= attn6;
					value += attn6*attn6*Extrapolate(xsb + 0, ysb + 1, zsb + 1, dx6, dy6, dz6);
				}
			}

			// First extra vertex
			double attnExt0 = 2 - dxExt0*dxExt0 - dyExt0*dyExt0 - dzExt0*dzExt0;
			if (attnExt0 > 0)
			{
				attnExt0 *= attnExt0;
				value += attnExt0*attnExt0*Extrapolate(xsvExt0, ysvExt0, zsvExt0, dxExt0, dyExt0, dzExt0);
			}

			// Second extra vertex
			double attnExt1 = 2 - dxExt1*dxExt1 - dyExt1*dyExt1 - dzExt1*dzExt1;
			if (attnExt1 > 0)
			{
				attnExt1 *= attnExt1;
				value += attnExt1*attnExt1*Extrapolate(xsvExt1, ysvExt1, zsvExt1, dxExt1, dyExt1, dzExt1);
			}

			return value/Normalization3D;
		}

		public double ComputeNoise(double x, double y, double z, double w)
		{
			// Place input coordinates on simplectic honeycomb.
			double stretchOffset = (x + y + z + w)*Stretch4D;
			double xs = x + stretchOffset;
			double ys = y + stretchOffset;
			double zs = z + stretchOffset;
			double ws = w + stretchOffset;

			// Floor to get simplectic honeycomb coordinates of rhombo-hypercube super-cell origin.
			int xsb = FastFloor(xs);
			int ysb = FastFloor(ys);
			int zsb = FastFloor(zs);
			int wsb = FastFloor(ws);

			// Skew out to get actual coordinates of stretched rhombo-hypercube origin. We'll need these later.
			double squishOffset = (xsb + ysb + zsb + wsb)*Squish4D;
			double xb = xsb + squishOffset;
			double yb = ysb + squishOffset;
			double zb = zsb + squishOffset;
			double wb = wsb + squishOffset;

			// Compute simplectic honeycomb coordinates relative to rhombo-hypercube origin.
			double xins = xs - xsb;
			double yins = ys - ysb;
			double zins = zs - zsb;
			double wins = ws - wsb;

			// Sum those together to get a value that determines which region we're in.
			double inSum = xins + yins + zins + wins;

			// Positions relative to origin point.
			double dx0 = x - xb;
			double dy0 = y - yb;
			double dz0 = z - zb;
			double dw0 = w - wb;

			// We'll be defining these inside the next block and using them afterwards.
			double dxExt0, dyExt0, dzExt0, dwExt0;
			double dxExt1, dyExt1, dzExt1, dwExt1;
			double dxExt2, dyExt2, dzExt2, dwExt2;
			int xsvExt0, ysvExt0, zsvExt0, wsvExt0;
			int xsvExt1, ysvExt1, zsvExt1, wsvExt1;
			int xsvExt2, ysvExt2, zsvExt2, wsvExt2;

			double value = 0;
			if (inSum <= 1)
			{
				// We're inside the pentachoron (4-Simplex) at (0,0,0,0)

				// Determine which two of (0,0,0,1), (0,0,1,0), (0,1,0,0), (1,0,0,0) are closest.
				sbyte aPoint = 0x01;
				double aScore = xins;
				sbyte bPoint = 0x02;
				double bScore = yins;
				if (aScore >= bScore && zins > bScore)
				{
					bScore = zins;
					bPoint = 0x04;
				}
				else if (aScore < bScore && zins > aScore)
				{
					aScore = zins;
					aPoint = 0x04;
				}
				if (aScore >= bScore && wins > bScore)
				{
					bScore = wins;
					bPoint = 0x08;
				}
				else if (aScore < bScore && wins > aScore)
				{
					aScore = wins;
					aPoint = 0x08;
				}

				// Now we determine the three lattice points not part of the pentachoron that may contribute.
				// This depends on the closest two pentachoron vertices, including (0,0,0,0)
				double uins = 1 - inSum;
				if (uins > aScore || uins > bScore)
				{
					// (0,0,0,0) is one of the closest two pentachoron vertices.
					sbyte c = (bScore > aScore ? bPoint : aPoint); // Our other closest vertex is the closest out of a and b.
					if ((c & 0x01) == 0)
					{
						xsvExt0 = xsb - 1;
						xsvExt1 = xsvExt2 = xsb;
						dxExt0 = dx0 + 1;
						dxExt1 = dxExt2 = dx0;
					}
					else
					{
						xsvExt0 = xsvExt1 = xsvExt2 = xsb + 1;
						dxExt0 = dxExt1 = dxExt2 = dx0 - 1;
					}

					if ((c & 0x02) == 0)
					{
						ysvExt0 = ysvExt1 = ysvExt2 = ysb;
						dyExt0 = dyExt1 = dyExt2 = dy0;
						if ((c & 0x01) == 0x01)
						{
							ysvExt0 -= 1;
							dyExt0 += 1;
						}
						else
						{
							ysvExt1 -= 1;
							dyExt1 += 1;
						}
					}
					else
					{
						ysvExt0 = ysvExt1 = ysvExt2 = ysb + 1;
						dyExt0 = dyExt1 = dyExt2 = dy0 - 1;
					}

					if ((c & 0x04) == 0)
					{
						zsvExt0 = zsvExt1 = zsvExt2 = zsb;
						dzExt0 = dzExt1 = dzExt2 = dz0;
						if ((c & 0x03) != 0)
						{
							if ((c & 0x03) == 0x03)
							{
								zsvExt0 -= 1;
								dzExt0 += 1;
							}
							else
							{
								zsvExt1 -= 1;
								dzExt1 += 1;
							}
						}
						else
						{
							zsvExt2 -= 1;
							dzExt2 += 1;
						}
					}
					else
					{
						zsvExt0 = zsvExt1 = zsvExt2 = zsb + 1;
						dzExt0 = dzExt1 = dzExt2 = dz0 - 1;
					}

					if ((c & 0x08) == 0)
					{
						wsvExt0 = wsvExt1 = wsb;
						wsvExt2 = wsb - 1;
						dwExt0 = dwExt1 = dw0;
						dwExt2 = dw0 + 1;
					}
					else
					{
						wsvExt0 = wsvExt1 = wsvExt2 = wsb + 1;
						dwExt0 = dwExt1 = dwExt2 = dw0 - 1;
					}
				}
				else
				{
					// (0,0,0,0) is not one of the closest two pentachoron vertices.
					var c = (sbyte) (aPoint | bPoint); // Our three extra vertices are determined by the closest two.

					if ((c & 0x01) == 0)
					{
						xsvExt0 = xsvExt2 = xsb;
						xsvExt1 = xsb - 1;
						dxExt0 = dx0 - 2*Squish4D;
						dxExt1 = dx0 + 1 - Squish4D;
						dxExt2 = dx0 - Squish4D;
					}
					else
					{
						xsvExt0 = xsvExt1 = xsvExt2 = xsb + 1;
						dxExt0 = dx0 - 1 - 2*Squish4D;
						dxExt1 = dxExt2 = dx0 - 1 - Squish4D;
					}

					if ((c & 0x02) == 0)
					{
						ysvExt0 = ysvExt1 = ysvExt2 = ysb;
						dyExt0 = dy0 - 2*Squish4D;
						dyExt1 = dyExt2 = dy0 - Squish4D;
						if ((c & 0x01) == 0x01)
						{
							ysvExt1 -= 1;
							dyExt1 += 1;
						}
						else
						{
							ysvExt2 -= 1;
							dyExt2 += 1;
						}
					}
					else
					{
						ysvExt0 = ysvExt1 = ysvExt2 = ysb + 1;
						dyExt0 = dy0 - 1 - 2*Squish4D;
						dyExt1 = dyExt2 = dy0 - 1 - Squish4D;
					}

					if ((c & 0x04) == 0)
					{
						zsvExt0 = zsvExt1 = zsvExt2 = zsb;
						dzExt0 = dz0 - 2*Squish4D;
						dzExt1 = dzExt2 = dz0 - Squish4D;
						if ((c & 0x03) == 0x03)
						{
							zsvExt1 -= 1;
							dzExt1 += 1;
						}
						else
						{
							zsvExt2 -= 1;
							dzExt2 += 1;
						}
					}
					else
					{
						zsvExt0 = zsvExt1 = zsvExt2 = zsb + 1;
						dzExt0 = dz0 - 1 - 2*Squish4D;
						dzExt1 = dzExt2 = dz0 - 1 - Squish4D;
					}

					if ((c & 0x08) == 0)
					{
						wsvExt0 = wsvExt1 = wsb;
						wsvExt2 = wsb - 1;
						dwExt0 = dw0 - 2*Squish4D;
						dwExt1 = dw0 - Squish4D;
						dwExt2 = dw0 + 1 - Squish4D;
					}
					else
					{
						wsvExt0 = wsvExt1 = wsvExt2 = wsb + 1;
						dwExt0 = dw0 - 1 - 2*Squish4D;
						dwExt1 = dwExt2 = dw0 - 1 - Squish4D;
					}
				}

				// Contribution (0,0,0,0)
				double attn0 = 2 - dx0*dx0 - dy0*dy0 - dz0*dz0 - dw0*dw0;
				if (attn0 > 0)
				{
					attn0 *= attn0;
					value += attn0*attn0*Extrapolate(xsb + 0, ysb + 0, zsb + 0, wsb + 0, dx0, dy0, dz0, dw0);
				}

				// Contribution (1,0,0,0)
				double dx1 = dx0 - 1 - Squish4D;
				double dy1 = dy0 - 0 - Squish4D;
				double dz1 = dz0 - 0 - Squish4D;
				double dw1 = dw0 - 0 - Squish4D;
				double attn1 = 2 - dx1*dx1 - dy1*dy1 - dz1*dz1 - dw1*dw1;
				if (attn1 > 0)
				{
					attn1 *= attn1;
					value += attn1*attn1*Extrapolate(xsb + 1, ysb + 0, zsb + 0, wsb + 0, dx1, dy1, dz1, dw1);
				}

				// Contribution (0,1,0,0)
				double dx2 = dx0 - 0 - Squish4D;
				double dy2 = dy0 - 1 - Squish4D;
				double dz2 = dz1;
				double dw2 = dw1;
				double attn2 = 2 - dx2*dx2 - dy2*dy2 - dz2*dz2 - dw2*dw2;
				if (attn2 > 0)
				{
					attn2 *= attn2;
					value += attn2*attn2*Extrapolate(xsb + 0, ysb + 1, zsb + 0, wsb + 0, dx2, dy2, dz2, dw2);
				}

				// Contribution (0,0,1,0)
				double dx3 = dx2;
				double dy3 = dy1;
				double dz3 = dz0 - 1 - Squish4D;
				double dw3 = dw1;
				double attn3 = 2 - dx3*dx3 - dy3*dy3 - dz3*dz3 - dw3*dw3;
				if (attn3 > 0)
				{
					attn3 *= attn3;
					value += attn3*attn3*Extrapolate(xsb + 0, ysb + 0, zsb + 1, wsb + 0, dx3, dy3, dz3, dw3);
				}

				// Contribution (0,0,0,1)
				double dx4 = dx2;
				double dy4 = dy1;
				double dz4 = dz1;
				double dw4 = dw0 - 1 - Squish4D;
				double attn4 = 2 - dx4*dx4 - dy4*dy4 - dz4*dz4 - dw4*dw4;
				if (attn4 > 0)
				{
					attn4 *= attn4;
					value += attn4*attn4*Extrapolate(xsb + 0, ysb + 0, zsb + 0, wsb + 1, dx4, dy4, dz4, dw4);
				}
			}
			else if (inSum >= 3)
			{
				// We're inside the pentachoron (4-Simplex) at (1,1,1,1)
				// Determine which two of (1,1,1,0), (1,1,0,1), (1,0,1,1), (0,1,1,1) are closest.
				sbyte aPoint = 0x0E;
				double aScore = xins;
				sbyte bPoint = 0x0D;
				double bScore = yins;
				if (aScore <= bScore && zins < bScore)
				{
					bScore = zins;
					bPoint = 0x0B;
				}
				else if (aScore > bScore && zins < aScore)
				{
					aScore = zins;
					aPoint = 0x0B;
				}
				if (aScore <= bScore && wins < bScore)
				{
					bScore = wins;
					bPoint = 0x07;
				}
				else if (aScore > bScore && wins < aScore)
				{
					aScore = wins;
					aPoint = 0x07;
				}

				// Now we determine the three lattice points not part of the pentachoron that may contribute.
				// This depends on the closest two pentachoron vertices, including (0,0,0,0)
				double uins = 4 - inSum;
				if (uins < aScore || uins < bScore)
				{
					// (1,1,1,1) is one of the closest two pentachoron vertices.
					sbyte c = (bScore < aScore ? bPoint : aPoint); // Our other closest vertex is the closest out of a and b.

					if ((c & 0x01) != 0)
					{
						xsvExt0 = xsb + 2;
						xsvExt1 = xsvExt2 = xsb + 1;
						dxExt0 = dx0 - 2 - 4*Squish4D;
						dxExt1 = dxExt2 = dx0 - 1 - 4*Squish4D;
					}
					else
					{
						xsvExt0 = xsvExt1 = xsvExt2 = xsb;
						dxExt0 = dxExt1 = dxExt2 = dx0 - 4*Squish4D;
					}

					if ((c & 0x02) != 0)
					{
						ysvExt0 = ysvExt1 = ysvExt2 = ysb + 1;
						dyExt0 = dyExt1 = dyExt2 = dy0 - 1 - 4*Squish4D;
						if ((c & 0x01) != 0)
						{
							ysvExt1 += 1;
							dyExt1 -= 1;
						}
						else
						{
							ysvExt0 += 1;
							dyExt0 -= 1;
						}
					}
					else
					{
						ysvExt0 = ysvExt1 = ysvExt2 = ysb;
						dyExt0 = dyExt1 = dyExt2 = dy0 - 4*Squish4D;
					}

					if ((c & 0x04) != 0)
					{
						zsvExt0 = zsvExt1 = zsvExt2 = zsb + 1;
						dzExt0 = dzExt1 = dzExt2 = dz0 - 1 - 4*Squish4D;
						if ((c & 0x03) != 0x03)
						{
							if ((c & 0x03) == 0)
							{
								zsvExt0 += 1;
								dzExt0 -= 1;
							}
							else
							{
								zsvExt1 += 1;
								dzExt1 -= 1;
							}
						}
						else
						{
							zsvExt2 += 1;
							dzExt2 -= 1;
						}
					}
					else
					{
						zsvExt0 = zsvExt1 = zsvExt2 = zsb;
						dzExt0 = dzExt1 = dzExt2 = dz0 - 4*Squish4D;
					}

					if ((c & 0x08) != 0)
					{
						wsvExt0 = wsvExt1 = wsb + 1;
						wsvExt2 = wsb + 2;
						dwExt0 = dwExt1 = dw0 - 1 - 4*Squish4D;
						dwExt2 = dw0 - 2 - 4*Squish4D;
					}
					else
					{
						wsvExt0 = wsvExt1 = wsvExt2 = wsb;
						dwExt0 = dwExt1 = dwExt2 = dw0 - 4*Squish4D;
					}
				}
				else
				{
					// (1,1,1,1) is not one of the closest two pentachoron vertices.
					var c = (sbyte) (aPoint & bPoint); // Our three extra vertices are determined by the closest two.

					if ((c & 0x01) != 0)
					{
						xsvExt0 = xsvExt2 = xsb + 1;
						xsvExt1 = xsb + 2;
						dxExt0 = dx0 - 1 - 2*Squish4D;
						dxExt1 = dx0 - 2 - 3*Squish4D;
						dxExt2 = dx0 - 1 - 3*Squish4D;
					}
					else
					{
						xsvExt0 = xsvExt1 = xsvExt2 = xsb;
						dxExt0 = dx0 - 2*Squish4D;
						dxExt1 = dxExt2 = dx0 - 3*Squish4D;
					}

					if ((c & 0x02) != 0)
					{
						ysvExt0 = ysvExt1 = ysvExt2 = ysb + 1;
						dyExt0 = dy0 - 1 - 2*Squish4D;
						dyExt1 = dyExt2 = dy0 - 1 - 3*Squish4D;
						if ((c & 0x01) != 0)
						{
							ysvExt2 += 1;
							dyExt2 -= 1;
						}
						else
						{
							ysvExt1 += 1;
							dyExt1 -= 1;
						}
					}
					else
					{
						ysvExt0 = ysvExt1 = ysvExt2 = ysb;
						dyExt0 = dy0 - 2*Squish4D;
						dyExt1 = dyExt2 = dy0 - 3*Squish4D;
					}

					if ((c & 0x04) != 0)
					{
						zsvExt0 = zsvExt1 = zsvExt2 = zsb + 1;
						dzExt0 = dz0 - 1 - 2*Squish4D;
						dzExt1 = dzExt2 = dz0 - 1 - 3*Squish4D;
						if ((c & 0x03) != 0)
						{
							zsvExt2 += 1;
							dzExt2 -= 1;
						}
						else
						{
							zsvExt1 += 1;
							dzExt1 -= 1;
						}
					}
					else
					{
						zsvExt0 = zsvExt1 = zsvExt2 = zsb;
						dzExt0 = dz0 - 2*Squish4D;
						dzExt1 = dzExt2 = dz0 - 3*Squish4D;
					}

					if ((c & 0x08) != 0)
					{
						wsvExt0 = wsvExt1 = wsb + 1;
						wsvExt2 = wsb + 2;
						dwExt0 = dw0 - 1 - 2*Squish4D;
						dwExt1 = dw0 - 1 - 3*Squish4D;
						dwExt2 = dw0 - 2 - 3*Squish4D;
					}
					else
					{
						wsvExt0 = wsvExt1 = wsvExt2 = wsb;
						dwExt0 = dw0 - 2*Squish4D;
						dwExt1 = dwExt2 = dw0 - 3*Squish4D;
					}
				}

				// Contribution (1,1,1,0)
				double dx4 = dx0 - 1 - 3*Squish4D;
				double dy4 = dy0 - 1 - 3*Squish4D;
				double dz4 = dz0 - 1 - 3*Squish4D;
				double dw4 = dw0 - 3*Squish4D;
				double attn4 = 2 - dx4*dx4 - dy4*dy4 - dz4*dz4 - dw4*dw4;
				if (attn4 > 0)
				{
					attn4 *= attn4;
					value += attn4*attn4*Extrapolate(xsb + 1, ysb + 1, zsb + 1, wsb + 0, dx4, dy4, dz4, dw4);
				}

				// Contribution (1,1,0,1)
				double dx3 = dx4;
				double dy3 = dy4;
				double dz3 = dz0 - 3*Squish4D;
				double dw3 = dw0 - 1 - 3*Squish4D;
				double attn3 = 2 - dx3*dx3 - dy3*dy3 - dz3*dz3 - dw3*dw3;
				if (attn3 > 0)
				{
					attn3 *= attn3;
					value += attn3*attn3*Extrapolate(xsb + 1, ysb + 1, zsb + 0, wsb + 1, dx3, dy3, dz3, dw3);
				}

				// Contribution (1,0,1,1)
				double dx2 = dx4;
				double dy2 = dy0 - 3*Squish4D;
				double dz2 = dz4;
				double dw2 = dw3;
				double attn2 = 2 - dx2*dx2 - dy2*dy2 - dz2*dz2 - dw2*dw2;
				if (attn2 > 0)
				{
					attn2 *= attn2;
					value += attn2*attn2*Extrapolate(xsb + 1, ysb + 0, zsb + 1, wsb + 1, dx2, dy2, dz2, dw2);
				}

				// Contribution (0,1,1,1)
				double dx1 = dx0 - 3*Squish4D;
				double dz1 = dz4;
				double dy1 = dy4;
				double dw1 = dw3;
				double attn1 = 2 - dx1*dx1 - dy1*dy1 - dz1*dz1 - dw1*dw1;
				if (attn1 > 0)
				{
					attn1 *= attn1;
					value += attn1*attn1*Extrapolate(xsb + 0, ysb + 1, zsb + 1, wsb + 1, dx1, dy1, dz1, dw1);
				}

				// Contribution (1,1,1,1)
				dx0 = dx0 - 1 - 4*Squish4D;
				dy0 = dy0 - 1 - 4*Squish4D;
				dz0 = dz0 - 1 - 4*Squish4D;
				dw0 = dw0 - 1 - 4*Squish4D;
				double attn0 = 2 - dx0*dx0 - dy0*dy0 - dz0*dz0 - dw0*dw0;
				if (attn0 > 0)
				{
					attn0 *= attn0;
					value += attn0*attn0*Extrapolate(xsb + 1, ysb + 1, zsb + 1, wsb + 1, dx0, dy0, dz0, dw0);
				}
			}
			else if (inSum <= 2)
			{
				// We're inside the first dispentachoron (Rectified 4-Simplex)
				double aScore;
				sbyte aPoint;
				bool aIsBiggerSide = true;
				double bScore;
				sbyte bPoint;
				bool bIsBiggerSide = true;

				// Decide between (1,1,0,0) and (0,0,1,1)
				if (xins + yins > zins + wins)
				{
					aScore = xins + yins;
					aPoint = 0x03;
				}
				else
				{
					aScore = zins + wins;
					aPoint = 0x0C;
				}

				// Decide between (1,0,1,0) and (0,1,0,1)
				if (xins + zins > yins + wins)
				{
					bScore = xins + zins;
					bPoint = 0x05;
				}
				else
				{
					bScore = yins + wins;
					bPoint = 0x0A;
				}

				// Closer between (1,0,0,1) and (0,1,1,0) will replace the further of a and b, if closer.
				if (xins + wins > yins + zins)
				{
					double score = xins + wins;
					if (aScore >= bScore && score > bScore)
					{
						bScore = score;
						bPoint = 0x09;
					}
					else if (aScore < bScore && score > aScore)
					{
						aScore = score;
						aPoint = 0x09;
					}
				}
				else
				{
					double score = yins + zins;
					if (aScore >= bScore && score > bScore)
					{
						bScore = score;
						bPoint = 0x06;
					}
					else if (aScore < bScore && score > aScore)
					{
						aScore = score;
						aPoint = 0x06;
					}
				}

				// Decide if (1,0,0,0) is closer.
				double p1 = 2 - inSum + xins;
				if (aScore >= bScore && p1 > bScore)
				{
					bScore = p1;
					bPoint = 0x01;
					bIsBiggerSide = false;
				}
				else if (aScore < bScore && p1 > aScore)
				{
					aScore = p1;
					aPoint = 0x01;
					aIsBiggerSide = false;
				}

				// Decide if (0,1,0,0) is closer.
				double p2 = 2 - inSum + yins;
				if (aScore >= bScore && p2 > bScore)
				{
					bScore = p2;
					bPoint = 0x02;
					bIsBiggerSide = false;
				}
				else if (aScore < bScore && p2 > aScore)
				{
					aScore = p2;
					aPoint = 0x02;
					aIsBiggerSide = false;
				}

				// Decide if (0,0,1,0) is closer.
				double p3 = 2 - inSum + zins;
				if (aScore >= bScore && p3 > bScore)
				{
					bScore = p3;
					bPoint = 0x04;
					bIsBiggerSide = false;
				}
				else if (aScore < bScore && p3 > aScore)
				{
					aScore = p3;
					aPoint = 0x04;
					aIsBiggerSide = false;
				}

				// Decide if (0,0,0,1) is closer.
				double p4 = 2 - inSum + wins;
				if (aScore >= bScore && p4 > bScore)
				{
					// bScore = p4;
					bPoint = 0x08;
					bIsBiggerSide = false;
				}
				else if (aScore < bScore && p4 > aScore)
				{
					// aScore = p4;
					aPoint = 0x08;
					aIsBiggerSide = false;
				}

				// Where each of the two closest points are determines how the extra three vertices are calculated.
				if (aIsBiggerSide == bIsBiggerSide)
				{
					if (aIsBiggerSide)
					{
						// Both closest points on the bigger side
						var c1 = (sbyte) (aPoint | bPoint);
						var c2 = (sbyte) (aPoint & bPoint);
						if ((c1 & 0x01) == 0)
						{
							xsvExt0 = xsb;
							xsvExt1 = xsb - 1;
							dxExt0 = dx0 - 3*Squish4D;
							dxExt1 = dx0 + 1 - 2*Squish4D;
						}
						else
						{
							xsvExt0 = xsvExt1 = xsb + 1;
							dxExt0 = dx0 - 1 - 3*Squish4D;
							dxExt1 = dx0 - 1 - 2*Squish4D;
						}

						if ((c1 & 0x02) == 0)
						{
							ysvExt0 = ysb;
							ysvExt1 = ysb - 1;
							dyExt0 = dy0 - 3*Squish4D;
							dyExt1 = dy0 + 1 - 2*Squish4D;
						}
						else
						{
							ysvExt0 = ysvExt1 = ysb + 1;
							dyExt0 = dy0 - 1 - 3*Squish4D;
							dyExt1 = dy0 - 1 - 2*Squish4D;
						}

						if ((c1 & 0x04) == 0)
						{
							zsvExt0 = zsb;
							zsvExt1 = zsb - 1;
							dzExt0 = dz0 - 3*Squish4D;
							dzExt1 = dz0 + 1 - 2*Squish4D;
						}
						else
						{
							zsvExt0 = zsvExt1 = zsb + 1;
							dzExt0 = dz0 - 1 - 3*Squish4D;
							dzExt1 = dz0 - 1 - 2*Squish4D;
						}

						if ((c1 & 0x08) == 0)
						{
							wsvExt0 = wsb;
							wsvExt1 = wsb - 1;
							dwExt0 = dw0 - 3*Squish4D;
							dwExt1 = dw0 + 1 - 2*Squish4D;
						}
						else
						{
							wsvExt0 = wsvExt1 = wsb + 1;
							dwExt0 = dw0 - 1 - 3*Squish4D;
							dwExt1 = dw0 - 1 - 2*Squish4D;
						}

						// One combination is a permutation of (0,0,0,2) based on c2
						xsvExt2 = xsb;
						ysvExt2 = ysb;
						zsvExt2 = zsb;
						wsvExt2 = wsb;
						dxExt2 = dx0 - 2*Squish4D;
						dyExt2 = dy0 - 2*Squish4D;
						dzExt2 = dz0 - 2*Squish4D;
						dwExt2 = dw0 - 2*Squish4D;
						if ((c2 & 0x01) != 0)
						{
							xsvExt2 += 2;
							dxExt2 -= 2;
						}
						else if ((c2 & 0x02) != 0)
						{
							ysvExt2 += 2;
							dyExt2 -= 2;
						}
						else if ((c2 & 0x04) != 0)
						{
							zsvExt2 += 2;
							dzExt2 -= 2;
						}
						else
						{
							wsvExt2 += 2;
							dwExt2 -= 2;
						}
					}
					else
					{
						// Both closest points on the smaller side
						// One of the two extra points is (0,0,0,0)
						xsvExt2 = xsb;
						ysvExt2 = ysb;
						zsvExt2 = zsb;
						wsvExt2 = wsb;
						dxExt2 = dx0;
						dyExt2 = dy0;
						dzExt2 = dz0;
						dwExt2 = dw0;

						// Other two points are based on the omitted axes.
						var c = (sbyte) (aPoint | bPoint);

						if ((c & 0x01) == 0)
						{
							xsvExt0 = xsb - 1;
							xsvExt1 = xsb;
							dxExt0 = dx0 + 1 - Squish4D;
							dxExt1 = dx0 - Squish4D;
						}
						else
						{
							xsvExt0 = xsvExt1 = xsb + 1;
							dxExt0 = dxExt1 = dx0 - 1 - Squish4D;
						}

						if ((c & 0x02) == 0)
						{
							ysvExt0 = ysvExt1 = ysb;
							dyExt0 = dyExt1 = dy0 - Squish4D;
							if ((c & 0x01) == 0x01)
							{
								ysvExt0 -= 1;
								dyExt0 += 1;
							}
							else
							{
								ysvExt1 -= 1;
								dyExt1 += 1;
							}
						}
						else
						{
							ysvExt0 = ysvExt1 = ysb + 1;
							dyExt0 = dyExt1 = dy0 - 1 - Squish4D;
						}

						if ((c & 0x04) == 0)
						{
							zsvExt0 = zsvExt1 = zsb;
							dzExt0 = dzExt1 = dz0 - Squish4D;
							if ((c & 0x03) == 0x03)
							{
								zsvExt0 -= 1;
								dzExt0 += 1;
							}
							else
							{
								zsvExt1 -= 1;
								dzExt1 += 1;
							}
						}
						else
						{
							zsvExt0 = zsvExt1 = zsb + 1;
							dzExt0 = dzExt1 = dz0 - 1 - Squish4D;
						}

						if ((c & 0x08) == 0)
						{
							wsvExt0 = wsb;
							wsvExt1 = wsb - 1;
							dwExt0 = dw0 - Squish4D;
							dwExt1 = dw0 + 1 - Squish4D;
						}
						else
						{
							wsvExt0 = wsvExt1 = wsb + 1;
							dwExt0 = dwExt1 = dw0 - 1 - Squish4D;
						}
					}
				}
				else
				{
					// One point on each "side"
					sbyte c1, c2;
					if (aIsBiggerSide)
					{
						c1 = aPoint;
						c2 = bPoint;
					}
					else
					{
						c1 = bPoint;
						c2 = aPoint;
					}

					// Two contributions are the bigger-sided point with each 0 replaced with -1.
					if ((c1 & 0x01) == 0)
					{
						xsvExt0 = xsb - 1;
						xsvExt1 = xsb;
						dxExt0 = dx0 + 1 - Squish4D;
						dxExt1 = dx0 - Squish4D;
					}
					else
					{
						xsvExt0 = xsvExt1 = xsb + 1;
						dxExt0 = dxExt1 = dx0 - 1 - Squish4D;
					}

					if ((c1 & 0x02) == 0)
					{
						ysvExt0 = ysvExt1 = ysb;
						dyExt0 = dyExt1 = dy0 - Squish4D;
						if ((c1 & 0x01) == 0x01)
						{
							ysvExt0 -= 1;
							dyExt0 += 1;
						}
						else
						{
							ysvExt1 -= 1;
							dyExt1 += 1;
						}
					}
					else
					{
						ysvExt0 = ysvExt1 = ysb + 1;
						dyExt0 = dyExt1 = dy0 - 1 - Squish4D;
					}

					if ((c1 & 0x04) == 0)
					{
						zsvExt0 = zsvExt1 = zsb;
						dzExt0 = dzExt1 = dz0 - Squish4D;
						if ((c1 & 0x03) == 0x03)
						{
							zsvExt0 -= 1;
							dzExt0 += 1;
						}
						else
						{
							zsvExt1 -= 1;
							dzExt1 += 1;
						}
					}
					else
					{
						zsvExt0 = zsvExt1 = zsb + 1;
						dzExt0 = dzExt1 = dz0 - 1 - Squish4D;
					}

					if ((c1 & 0x08) == 0)
					{
						wsvExt0 = wsb;
						wsvExt1 = wsb - 1;
						dwExt0 = dw0 - Squish4D;
						dwExt1 = dw0 + 1 - Squish4D;
					}
					else
					{
						wsvExt0 = wsvExt1 = wsb + 1;
						dwExt0 = dwExt1 = dw0 - 1 - Squish4D;
					}

					// One contribution is a permutation of (0,0,0,2) based on the smaller-sided point
					xsvExt2 = xsb;
					ysvExt2 = ysb;
					zsvExt2 = zsb;
					wsvExt2 = wsb;
					dxExt2 = dx0 - 2*Squish4D;
					dyExt2 = dy0 - 2*Squish4D;
					dzExt2 = dz0 - 2*Squish4D;
					dwExt2 = dw0 - 2*Squish4D;
					if ((c2 & 0x01) != 0)
					{
						xsvExt2 += 2;
						dxExt2 -= 2;
					}
					else if ((c2 & 0x02) != 0)
					{
						ysvExt2 += 2;
						dyExt2 -= 2;
					}
					else if ((c2 & 0x04) != 0)
					{
						zsvExt2 += 2;
						dzExt2 -= 2;
					}
					else
					{
						wsvExt2 += 2;
						dwExt2 -= 2;
					}
				}

				// Contribution (1,0,0,0)
				double dx1 = dx0 - 1 - Squish4D;
				double dy1 = dy0 - 0 - Squish4D;
				double dz1 = dz0 - 0 - Squish4D;
				double dw1 = dw0 - 0 - Squish4D;
				double attn1 = 2 - dx1*dx1 - dy1*dy1 - dz1*dz1 - dw1*dw1;
				if (attn1 > 0)
				{
					attn1 *= attn1;
					value += attn1*attn1*Extrapolate(xsb + 1, ysb + 0, zsb + 0, wsb + 0, dx1, dy1, dz1, dw1);
				}

				// Contribution (0,1,0,0)
				double dx2 = dx0 - 0 - Squish4D;
				double dy2 = dy0 - 1 - Squish4D;
				double dz2 = dz1;
				double dw2 = dw1;
				double attn2 = 2 - dx2*dx2 - dy2*dy2 - dz2*dz2 - dw2*dw2;
				if (attn2 > 0)
				{
					attn2 *= attn2;
					value += attn2*attn2*Extrapolate(xsb + 0, ysb + 1, zsb + 0, wsb + 0, dx2, dy2, dz2, dw2);
				}

				// Contribution (0,0,1,0)
				double dx3 = dx2;
				double dy3 = dy1;
				double dz3 = dz0 - 1 - Squish4D;
				double dw3 = dw1;
				double attn3 = 2 - dx3*dx3 - dy3*dy3 - dz3*dz3 - dw3*dw3;
				if (attn3 > 0)
				{
					attn3 *= attn3;
					value += attn3*attn3*Extrapolate(xsb + 0, ysb + 0, zsb + 1, wsb + 0, dx3, dy3, dz3, dw3);
				}

				// Contribution (0,0,0,1)
				double dx4 = dx2;
				double dy4 = dy1;
				double dz4 = dz1;
				double dw4 = dw0 - 1 - Squish4D;
				double attn4 = 2 - dx4*dx4 - dy4*dy4 - dz4*dz4 - dw4*dw4;
				if (attn4 > 0)
				{
					attn4 *= attn4;
					value += attn4*attn4*Extrapolate(xsb + 0, ysb + 0, zsb + 0, wsb + 1, dx4, dy4, dz4, dw4);
				}

				// Contribution (1,1,0,0)
				double dx5 = dx0 - 1 - 2*Squish4D;
				double dy5 = dy0 - 1 - 2*Squish4D;
				double dz5 = dz0 - 0 - 2*Squish4D;
				double dw5 = dw0 - 0 - 2*Squish4D;
				double attn5 = 2 - dx5*dx5 - dy5*dy5 - dz5*dz5 - dw5*dw5;
				if (attn5 > 0)
				{
					attn5 *= attn5;
					value += attn5*attn5*Extrapolate(xsb + 1, ysb + 1, zsb + 0, wsb + 0, dx5, dy5, dz5, dw5);
				}

				// Contribution (1,0,1,0)
				double dx6 = dx0 - 1 - 2*Squish4D;
				double dy6 = dy0 - 0 - 2*Squish4D;
				double dz6 = dz0 - 1 - 2*Squish4D;
				double dw6 = dw0 - 0 - 2*Squish4D;
				double attn6 = 2 - dx6*dx6 - dy6*dy6 - dz6*dz6 - dw6*dw6;
				if (attn6 > 0)
				{
					attn6 *= attn6;
					value += attn6*attn6*Extrapolate(xsb + 1, ysb + 0, zsb + 1, wsb + 0, dx6, dy6, dz6, dw6);
				}

				// Contribution (1,0,0,1)
				double dx7 = dx0 - 1 - 2*Squish4D;
				double dy7 = dy0 - 0 - 2*Squish4D;
				double dz7 = dz0 - 0 - 2*Squish4D;
				double dw7 = dw0 - 1 - 2*Squish4D;
				double attn7 = 2 - dx7*dx7 - dy7*dy7 - dz7*dz7 - dw7*dw7;
				if (attn7 > 0)
				{
					attn7 *= attn7;
					value += attn7*attn7*Extrapolate(xsb + 1, ysb + 0, zsb + 0, wsb + 1, dx7, dy7, dz7, dw7);
				}

				// Contribution (0,1,1,0)
				double dx8 = dx0 - 0 - 2*Squish4D;
				double dy8 = dy0 - 1 - 2*Squish4D;
				double dz8 = dz0 - 1 - 2*Squish4D;
				double dw8 = dw0 - 0 - 2*Squish4D;
				double attn8 = 2 - dx8*dx8 - dy8*dy8 - dz8*dz8 - dw8*dw8;
				if (attn8 > 0)
				{
					attn8 *= attn8;
					value += attn8*attn8*Extrapolate(xsb + 0, ysb + 1, zsb + 1, wsb + 0, dx8, dy8, dz8, dw8);
				}

				// Contribution (0,1,0,1)
				double dx9 = dx0 - 0 - 2*Squish4D;
				double dy9 = dy0 - 1 - 2*Squish4D;
				double dz9 = dz0 - 0 - 2*Squish4D;
				double dw9 = dw0 - 1 - 2*Squish4D;
				double attn9 = 2 - dx9*dx9 - dy9*dy9 - dz9*dz9 - dw9*dw9;
				if (attn9 > 0)
				{
					attn9 *= attn9;
					value += attn9*attn9*Extrapolate(xsb + 0, ysb + 1, zsb + 0, wsb + 1, dx9, dy9, dz9, dw9);
				}

				// Contribution (0,0,1,1)
				double dx10 = dx0 - 0 - 2*Squish4D;
				double dy10 = dy0 - 0 - 2*Squish4D;
				double dz10 = dz0 - 1 - 2*Squish4D;
				double dw10 = dw0 - 1 - 2*Squish4D;
				double attn10 = 2 - dx10*dx10 - dy10*dy10 - dz10*dz10 - dw10*dw10;
				if (attn10 > 0)
				{
					attn10 *= attn10;
					value += attn10*attn10*Extrapolate(xsb + 0, ysb + 0, zsb + 1, wsb + 1, dx10, dy10, dz10, dw10);
				}
			}
			else
			{
				// We're inside the second dispentachoron (Rectified 4-Simplex)
				double aScore;
				sbyte aPoint;
				bool aIsBiggerSide = true;
				double bScore;
				sbyte bPoint;
				bool bIsBiggerSide = true;

				// Decide between (0,0,1,1) and (1,1,0,0)
				if (xins + yins < zins + wins)
				{
					aScore = xins + yins;
					aPoint = 0x0C;
				}
				else
				{
					aScore = zins + wins;
					aPoint = 0x03;
				}

				// Decide between (0,1,0,1) and (1,0,1,0)
				if (xins + zins < yins + wins)
				{
					bScore = xins + zins;
					bPoint = 0x0A;
				}
				else
				{
					bScore = yins + wins;
					bPoint = 0x05;
				}

				// Closer between (0,1,1,0) and (1,0,0,1) will replace the further of a and b, if closer.
				if (xins + wins < yins + zins)
				{
					double score = xins + wins;
					if (aScore <= bScore && score < bScore)
					{
						bScore = score;
						bPoint = 0x06;
					}
					else if (aScore > bScore && score < aScore)
					{
						aScore = score;
						aPoint = 0x06;
					}
				}
				else
				{
					double score = yins + zins;
					if (aScore <= bScore && score < bScore)
					{
						bScore = score;
						bPoint = 0x09;
					}
					else if (aScore > bScore && score < aScore)
					{
						aScore = score;
						aPoint = 0x09;
					}
				}

				// Decide if (0,1,1,1) is closer.
				double p1 = 3 - inSum + xins;
				if (aScore <= bScore && p1 < bScore)
				{
					bScore = p1;
					bPoint = 0x0E;
					bIsBiggerSide = false;
				}
				else if (aScore > bScore && p1 < aScore)
				{
					aScore = p1;
					aPoint = 0x0E;
					aIsBiggerSide = false;
				}

				// Decide if (1,0,1,1) is closer.
				double p2 = 3 - inSum + yins;
				if (aScore <= bScore && p2 < bScore)
				{
					bScore = p2;
					bPoint = 0x0D;
					bIsBiggerSide = false;
				}
				else if (aScore > bScore && p2 < aScore)
				{
					aScore = p2;
					aPoint = 0x0D;
					aIsBiggerSide = false;
				}

				// Decide if (1,1,0,1) is closer.
				double p3 = 3 - inSum + zins;
				if (aScore <= bScore && p3 < bScore)
				{
					bScore = p3;
					bPoint = 0x0B;
					bIsBiggerSide = false;
				}
				else if (aScore > bScore && p3 < aScore)
				{
					aScore = p3;
					aPoint = 0x0B;
					aIsBiggerSide = false;
				}

				// Decide if (1,1,1,0) is closer.
				double p4 = 3 - inSum + wins;
				if (aScore <= bScore && p4 < bScore)
				{
					// bScore = p4;
					bPoint = 0x07;
					bIsBiggerSide = false;
				}
				else if (aScore > bScore && p4 < aScore)
				{
					// aScore = p4;
					aPoint = 0x07;
					aIsBiggerSide = false;
				}

				// Where each of the two closest points are determines how the extra three vertices are calculated.
				if (aIsBiggerSide == bIsBiggerSide)
				{
					if (aIsBiggerSide)
					{
						// Both closest points on the bigger side
						var c1 = (sbyte) (aPoint & bPoint);
						var c2 = (sbyte) (aPoint | bPoint);

						// Two contributions are permutations of (0,0,0,1) and (0,0,0,2) based on c1
						xsvExt0 = xsvExt1 = xsb;
						ysvExt0 = ysvExt1 = ysb;
						zsvExt0 = zsvExt1 = zsb;
						wsvExt0 = wsvExt1 = wsb;
						dxExt0 = dx0 - Squish4D;
						dyExt0 = dy0 - Squish4D;
						dzExt0 = dz0 - Squish4D;
						dwExt0 = dw0 - Squish4D;
						dxExt1 = dx0 - 2*Squish4D;
						dyExt1 = dy0 - 2*Squish4D;
						dzExt1 = dz0 - 2*Squish4D;
						dwExt1 = dw0 - 2*Squish4D;
						if ((c1 & 0x01) != 0)
						{
							xsvExt0 += 1;
							dxExt0 -= 1;
							xsvExt1 += 2;
							dxExt1 -= 2;
						}
						else if ((c1 & 0x02) != 0)
						{
							ysvExt0 += 1;
							dyExt0 -= 1;
							ysvExt1 += 2;
							dyExt1 -= 2;
						}
						else if ((c1 & 0x04) != 0)
						{
							zsvExt0 += 1;
							dzExt0 -= 1;
							zsvExt1 += 2;
							dzExt1 -= 2;
						}
						else
						{
							wsvExt0 += 1;
							dwExt0 -= 1;
							wsvExt1 += 2;
							dwExt1 -= 2;
						}

						// One contribution is a permutation of (1,1,1,-1) based on c2
						xsvExt2 = xsb + 1;
						ysvExt2 = ysb + 1;
						zsvExt2 = zsb + 1;
						wsvExt2 = wsb + 1;
						dxExt2 = dx0 - 1 - 2*Squish4D;
						dyExt2 = dy0 - 1 - 2*Squish4D;
						dzExt2 = dz0 - 1 - 2*Squish4D;
						dwExt2 = dw0 - 1 - 2*Squish4D;
						if ((c2 & 0x01) == 0)
						{
							xsvExt2 -= 2;
							dxExt2 += 2;
						}
						else if ((c2 & 0x02) == 0)
						{
							ysvExt2 -= 2;
							dyExt2 += 2;
						}
						else if ((c2 & 0x04) == 0)
						{
							zsvExt2 -= 2;
							dzExt2 += 2;
						}
						else
						{
							wsvExt2 -= 2;
							dwExt2 += 2;
						}
					}
					else
					{
						// Both closest points on the smaller side
						// One of the two extra points is (1,1,1,1)
						xsvExt2 = xsb + 1;
						ysvExt2 = ysb + 1;
						zsvExt2 = zsb + 1;
						wsvExt2 = wsb + 1;
						dxExt2 = dx0 - 1 - 4*Squish4D;
						dyExt2 = dy0 - 1 - 4*Squish4D;
						dzExt2 = dz0 - 1 - 4*Squish4D;
						dwExt2 = dw0 - 1 - 4*Squish4D;

						// Other two points are based on the shared axes.
						var c = (sbyte) (aPoint & bPoint);

						if ((c & 0x01) != 0)
						{
							xsvExt0 = xsb + 2;
							xsvExt1 = xsb + 1;
							dxExt0 = dx0 - 2 - 3*Squish4D;
							dxExt1 = dx0 - 1 - 3*Squish4D;
						}
						else
						{
							xsvExt0 = xsvExt1 = xsb;
							dxExt0 = dxExt1 = dx0 - 3*Squish4D;
						}

						if ((c & 0x02) != 0)
						{
							ysvExt0 = ysvExt1 = ysb + 1;
							dyExt0 = dyExt1 = dy0 - 1 - 3*Squish4D;
							if ((c & 0x01) == 0)
							{
								ysvExt0 += 1;
								dyExt0 -= 1;
							}
							else
							{
								ysvExt1 += 1;
								dyExt1 -= 1;
							}
						}
						else
						{
							ysvExt0 = ysvExt1 = ysb;
							dyExt0 = dyExt1 = dy0 - 3*Squish4D;
						}

						if ((c & 0x04) != 0)
						{
							zsvExt0 = zsvExt1 = zsb + 1;
							dzExt0 = dzExt1 = dz0 - 1 - 3*Squish4D;
							if ((c & 0x03) == 0)
							{
								zsvExt0 += 1;
								dzExt0 -= 1;
							}
							else
							{
								zsvExt1 += 1;
								dzExt1 -= 1;
							}
						}
						else
						{
							zsvExt0 = zsvExt1 = zsb;
							dzExt0 = dzExt1 = dz0 - 3*Squish4D;
						}

						if ((c & 0x08) != 0)
						{
							wsvExt0 = wsb + 1;
							wsvExt1 = wsb + 2;
							dwExt0 = dw0 - 1 - 3*Squish4D;
							dwExt1 = dw0 - 2 - 3*Squish4D;
						}
						else
						{
							wsvExt0 = wsvExt1 = wsb;
							dwExt0 = dwExt1 = dw0 - 3*Squish4D;
						}
					}
				}
				else
				{
					// One point on each "side"
					sbyte c1, c2;
					if (aIsBiggerSide)
					{
						c1 = aPoint;
						c2 = bPoint;
					}
					else
					{
						c1 = bPoint;
						c2 = aPoint;
					}

					// Two contributions are the bigger-sided point with each 1 replaced with 2.
					if ((c1 & 0x01) != 0)
					{
						xsvExt0 = xsb + 2;
						xsvExt1 = xsb + 1;
						dxExt0 = dx0 - 2 - 3*Squish4D;
						dxExt1 = dx0 - 1 - 3*Squish4D;
					}
					else
					{
						xsvExt0 = xsvExt1 = xsb;
						dxExt0 = dxExt1 = dx0 - 3*Squish4D;
					}

					if ((c1 & 0x02) != 0)
					{
						ysvExt0 = ysvExt1 = ysb + 1;
						dyExt0 = dyExt1 = dy0 - 1 - 3*Squish4D;
						if ((c1 & 0x01) == 0)
						{
							ysvExt0 += 1;
							dyExt0 -= 1;
						}
						else
						{
							ysvExt1 += 1;
							dyExt1 -= 1;
						}
					}
					else
					{
						ysvExt0 = ysvExt1 = ysb;
						dyExt0 = dyExt1 = dy0 - 3*Squish4D;
					}

					if ((c1 & 0x04) != 0)
					{
						zsvExt0 = zsvExt1 = zsb + 1;
						dzExt0 = dzExt1 = dz0 - 1 - 3*Squish4D;
						if ((c1 & 0x03) == 0)
						{
							zsvExt0 += 1;
							dzExt0 -= 1;
						}
						else
						{
							zsvExt1 += 1;
							dzExt1 -= 1;
						}
					}
					else
					{
						zsvExt0 = zsvExt1 = zsb;
						dzExt0 = dzExt1 = dz0 - 3*Squish4D;
					}

					if ((c1 & 0x08) != 0)
					{
						wsvExt0 = wsb + 1;
						wsvExt1 = wsb + 2;
						dwExt0 = dw0 - 1 - 3*Squish4D;
						dwExt1 = dw0 - 2 - 3*Squish4D;
					}
					else
					{
						wsvExt0 = wsvExt1 = wsb;
						dwExt0 = dwExt1 = dw0 - 3*Squish4D;
					}

					// One contribution is a permutation of (1,1,1,-1) based on the smaller-sided point
					xsvExt2 = xsb + 1;
					ysvExt2 = ysb + 1;
					zsvExt2 = zsb + 1;
					wsvExt2 = wsb + 1;
					dxExt2 = dx0 - 1 - 2*Squish4D;
					dyExt2 = dy0 - 1 - 2*Squish4D;
					dzExt2 = dz0 - 1 - 2*Squish4D;
					dwExt2 = dw0 - 1 - 2*Squish4D;
					if ((c2 & 0x01) == 0)
					{
						xsvExt2 -= 2;
						dxExt2 += 2;
					}
					else if ((c2 & 0x02) == 0)
					{
						ysvExt2 -= 2;
						dyExt2 += 2;
					}
					else if ((c2 & 0x04) == 0)
					{
						zsvExt2 -= 2;
						dzExt2 += 2;
					}
					else
					{
						wsvExt2 -= 2;
						dwExt2 += 2;
					}
				}

				// Contribution (1,1,1,0)
				double dx4 = dx0 - 1 - 3*Squish4D;
				double dy4 = dy0 - 1 - 3*Squish4D;
				double dz4 = dz0 - 1 - 3*Squish4D;
				double dw4 = dw0 - 3*Squish4D;
				double attn4 = 2 - dx4*dx4 - dy4*dy4 - dz4*dz4 - dw4*dw4;
				if (attn4 > 0)
				{
					attn4 *= attn4;
					value += attn4*attn4*Extrapolate(xsb + 1, ysb + 1, zsb + 1, wsb + 0, dx4, dy4, dz4, dw4);
				}

				// Contribution (1,1,0,1)
				double dx3 = dx4;
				double dy3 = dy4;
				double dz3 = dz0 - 3*Squish4D;
				double dw3 = dw0 - 1 - 3*Squish4D;
				double attn3 = 2 - dx3*dx3 - dy3*dy3 - dz3*dz3 - dw3*dw3;
				if (attn3 > 0)
				{
					attn3 *= attn3;
					value += attn3*attn3*Extrapolate(xsb + 1, ysb + 1, zsb + 0, wsb + 1, dx3, dy3, dz3, dw3);
				}

				// Contribution (1,0,1,1)
				double dx2 = dx4;
				double dy2 = dy0 - 3*Squish4D;
				double dz2 = dz4;
				double dw2 = dw3;
				double attn2 = 2 - dx2*dx2 - dy2*dy2 - dz2*dz2 - dw2*dw2;
				if (attn2 > 0)
				{
					attn2 *= attn2;
					value += attn2*attn2*Extrapolate(xsb + 1, ysb + 0, zsb + 1, wsb + 1, dx2, dy2, dz2, dw2);
				}

				// Contribution (0,1,1,1)
				double dx1 = dx0 - 3*Squish4D;
				double dz1 = dz4;
				double dy1 = dy4;
				double dw1 = dw3;
				double attn1 = 2 - dx1*dx1 - dy1*dy1 - dz1*dz1 - dw1*dw1;
				if (attn1 > 0)
				{
					attn1 *= attn1;
					value += attn1*attn1*Extrapolate(xsb + 0, ysb + 1, zsb + 1, wsb + 1, dx1, dy1, dz1, dw1);
				}

				// Contribution (1,1,0,0)
				double dx5 = dx0 - 1 - 2*Squish4D;
				double dy5 = dy0 - 1 - 2*Squish4D;
				double dz5 = dz0 - 0 - 2*Squish4D;
				double dw5 = dw0 - 0 - 2*Squish4D;
				double attn5 = 2 - dx5*dx5 - dy5*dy5 - dz5*dz5 - dw5*dw5;
				if (attn5 > 0)
				{
					attn5 *= attn5;
					value += attn5*attn5*Extrapolate(xsb + 1, ysb + 1, zsb + 0, wsb + 0, dx5, dy5, dz5, dw5);
				}

				// Contribution (1,0,1,0)
				double dx6 = dx0 - 1 - 2*Squish4D;
				double dy6 = dy0 - 0 - 2*Squish4D;
				double dz6 = dz0 - 1 - 2*Squish4D;
				double dw6 = dw0 - 0 - 2*Squish4D;
				double attn6 = 2 - dx6*dx6 - dy6*dy6 - dz6*dz6 - dw6*dw6;
				if (attn6 > 0)
				{
					attn6 *= attn6;
					value += attn6*attn6*Extrapolate(xsb + 1, ysb + 0, zsb + 1, wsb + 0, dx6, dy6, dz6, dw6);
				}

				// Contribution (1,0,0,1)
				double dx7 = dx0 - 1 - 2*Squish4D;
				double dy7 = dy0 - 0 - 2*Squish4D;
				double dz7 = dz0 - 0 - 2*Squish4D;
				double dw7 = dw0 - 1 - 2*Squish4D;
				double attn7 = 2 - dx7*dx7 - dy7*dy7 - dz7*dz7 - dw7*dw7;
				if (attn7 > 0)
				{
					attn7 *= attn7;
					value += attn7*attn7*Extrapolate(xsb + 1, ysb + 0, zsb + 0, wsb + 1, dx7, dy7, dz7, dw7);
				}

				// Contribution (0,1,1,0)
				double dx8 = dx0 - 0 - 2*Squish4D;
				double dy8 = dy0 - 1 - 2*Squish4D;
				double dz8 = dz0 - 1 - 2*Squish4D;
				double dw8 = dw0 - 0 - 2*Squish4D;
				double attn8 = 2 - dx8*dx8 - dy8*dy8 - dz8*dz8 - dw8*dw8;
				if (attn8 > 0)
				{
					attn8 *= attn8;
					value += attn8*attn8*Extrapolate(xsb + 0, ysb + 1, zsb + 1, wsb + 0, dx8, dy8, dz8, dw8);
				}

				// Contribution (0,1,0,1)
				double dx9 = dx0 - 0 - 2*Squish4D;
				double dy9 = dy0 - 1 - 2*Squish4D;
				double dz9 = dz0 - 0 - 2*Squish4D;
				double dw9 = dw0 - 1 - 2*Squish4D;
				double attn9 = 2 - dx9*dx9 - dy9*dy9 - dz9*dz9 - dw9*dw9;
				if (attn9 > 0)
				{
					attn9 *= attn9;
					value += attn9*attn9*Extrapolate(xsb + 0, ysb + 1, zsb + 0, wsb + 1, dx9, dy9, dz9, dw9);
				}

				// Contribution (0,0,1,1)
				double dx10 = dx0 - 0 - 2*Squish4D;
				double dy10 = dy0 - 0 - 2*Squish4D;
				double dz10 = dz0 - 1 - 2*Squish4D;
				double dw10 = dw0 - 1 - 2*Squish4D;
				double attn10 = 2 - dx10*dx10 - dy10*dy10 - dz10*dz10 - dw10*dw10;
				if (attn10 > 0)
				{
					attn10 *= attn10;
					value += attn10*attn10*Extrapolate(xsb + 0, ysb + 0, zsb + 1, wsb + 1, dx10, dy10, dz10, dw10);
				}
			}

			// First extra vertex
			double attnExt0 = 2 - dxExt0*dxExt0 - dyExt0*dyExt0 - dzExt0*dzExt0 - dwExt0*dwExt0;
			if (attnExt0 > 0)
			{
				attnExt0 *= attnExt0;
				value += attnExt0*attnExt0*Extrapolate(xsvExt0, ysvExt0, zsvExt0, wsvExt0, dxExt0, dyExt0, dzExt0, dwExt0);
			}

			// Second extra vertex
			double attnExt1 = 2 - dxExt1*dxExt1 - dyExt1*dyExt1 - dzExt1*dzExt1 - dwExt1*dwExt1;
			if (attnExt1 > 0)
			{
				attnExt1 *= attnExt1;
				value += attnExt1*attnExt1*Extrapolate(xsvExt1, ysvExt1, zsvExt1, wsvExt1, dxExt1, dyExt1, dzExt1, dwExt1);
			}

			// Third extra vertex
			double attnExt2 = 2 - dxExt2*dxExt2 - dyExt2*dyExt2 - dzExt2*dzExt2 - dwExt2*dwExt2;
			if (attnExt2 > 0)
			{
				attnExt2 *= attnExt2;
				value += attnExt2*attnExt2*Extrapolate(xsvExt2, ysvExt2, zsvExt2, wsvExt2, dxExt2, dyExt2, dzExt2, dwExt2);
			}

			return value/Normalization4D;
		}

		private double Extrapolate(int xsb, int ysb, double dx, double dy)
		{
			int index = _permutation[(_permutation[xsb & 0xFF] + ysb) & 0xFF] & 0x0E;
			return Gradients2D[index]*dx
			       + Gradients2D[index + 1]*dy;
		}

		private double Extrapolate(int xsb, int ysb, int zsb, double dx, double dy, double dz)
		{
			int index = _permutationGradIndex3D[(_permutation[(_permutation[xsb & 0xFF] + ysb) & 0xFF] + zsb) & 0xFF];
			return Gradients3D[index]*dx
			       + Gradients3D[index + 1]*dy
			       + Gradients3D[index + 2]*dz;
		}

		private double Extrapolate(int xsb, int ysb, int zsb, int wsb, double dx, double dy, double dz, double dw)
		{
			int index = _permutation[(_permutation[(_permutation[(_permutation[xsb & 0xFF]
			                                                      + ysb) & 0xFF] + zsb) & 0xFF] + wsb) & 0xFF] & 0xFC;
			return Gradients4D[index]*dx
			       + Gradients4D[index + 1]*dy
			       + Gradients4D[index + 2]*dz
			       + Gradients4D[index + 3]*dw;
		}

		private static int FastFloor(double x)
		{
			var xi = (int) x;
			return x < xi ? xi - 1 : xi;
		}
	}
}