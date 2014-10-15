namespace BitFn.NoiseLibrary.Algorithms
{
	/// <summary>
	///     A port of Kurt Spencer's SimplexValueNoise.java. Computes coherent noise based on value-multiplied kernels on a
	///     triangular grid.
	/// </summary>
	/// <remarks>
	///     Original Java implementation released under public domain:
	///     https://gist.github.com/KdotJPG/9bbab7d3655b82811b24/22f0b3424c44829846390235da258f5ea705e71c
	/// </remarks>
	public class SimplexValueNoise : INoiseService2D
	{
		/// <summary>
		///     The squish constant to use for two dimensional noise.
		/// </summary>
		/// <remarks>
		///     <code>(sqrt(2+1)-1)/2 = 0.366025403784438646763723170752936183471402626905190314027903</code>
		/// </remarks>
		private const double Squish = 0.366025403784438646763723170752936183471402626905190314027903;

		/// <summary>
		///     The strech constant to use for two dimensional noise.
		/// </summary>
		/// <remarks>
		///     <code>(1/sqrt(2+1)-1)/2 = -0.21132486540518711774542560974902127217619912436493656199069</code>
		/// </remarks>
		private const double Stretch = -0.21132486540518711774542560974902127217619912436493656199069;

		/// <summary>
		///     The normalization scalar to use for two dimensional noise.
		/// </summary>
		private const double Normalization = 2174;

		/// <summary>
		///     The default seed to use for the parameterless constructor.
		/// </summary>
		private const long DefaultSeed = 0;

		/// <summary>
		///     The radius of the two dimensional kernal.
		/// </summary>
		private const double KernelRadius = 1.732050807568877293527446341505872366942805253810380628055806;

		private const int PermutationLength = 256;
		private const long PermutationConstant1 = 6364136223846793005L;
		private const long PermutationConstant2 = 1442695040888963407L;

		private readonly byte[] _permutation;

		public SimplexValueNoise()
			: this(DefaultSeed)
		{
		}

		public SimplexValueNoise(byte[] perm)
		{
			_permutation = perm;
		}

		/// <summary>
		///     Initializes the class using a permutation array generated from a 64-bit seed. Generates a proper permutation (i.e.
		///     doesn't merely perform N successive pair swaps on a base array) Uses a simple 64-bit LCG.
		/// </summary>
		/// <param name="seed"></param>
		public SimplexValueNoise(long seed)
		{
			_permutation = new byte[PermutationLength];
			var source = new byte[PermutationLength];
			for (short i = 0; i < PermutationLength; i++)
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
				source[r] = source[i];
			}
		}

		public double ComputeNoise(double x, double y)
		{
			// Place input coordinates on triangular grid.
			double squishOffset = (x + y)*Squish;
			double xs = x + squishOffset;
			double ys = y + squishOffset;

			// Floor to get base coordinate of containing square/rhombus.
			int xsb = FastFloor(xs);
			int ysb = FastFloor(ys);

			// Skew out to get actual coordinates of rhombus origin. We'll need these later.
			double stretchOffset = (xsb + ysb)*Stretch;
			double xb = xsb + stretchOffset;
			double yb = ysb + stretchOffset;

			// Positions relative to origin point.
			double dx = x - xb;
			double dy = y - yb;

			// Compute grid coordinates relative to rhombus origin.
			double xins = xs - xsb;
			double yins = ys - ysb;

			double value;
			if (xins > yins)
			{
				// We're inside the x>y triangle of the rhombus

				// Get our 12 surrounding vertex values
				short yp = _permutation[(ysb - 1) & 0xFF];
				var h1 = (sbyte) _permutation[(yp + xsb - 1) & 0xFF]; // (-1,-1)
				var h2 = (sbyte) _permutation[(yp + xsb + 0) & 0xFF]; // ( 0,-1)
				var h3 = (sbyte) _permutation[(yp + xsb + 1) & 0xFF]; // ( 1,-1)
				yp = _permutation[(ysb + 0) & 0xFF];
				var h4 = (sbyte) _permutation[(yp + xsb - 1) & 0xFF]; // (-1, 0)
				var h5 = (sbyte) _permutation[(yp + xsb + 0) & 0xFF]; // ( 0, 0)
				var h6 = (sbyte) _permutation[(yp + xsb + 1) & 0xFF]; // ( 1, 0)
				var h7 = (sbyte) _permutation[(yp + xsb + 2) & 0xFF]; // ( 2, 0)
				yp = _permutation[(ysb + 1) & 0xFF];
				var h8 = (sbyte) _permutation[(yp + xsb + 0) & 0xFF]; // ( 0, 1)
				var h9 = (sbyte) _permutation[(yp + xsb + 1) & 0xFF]; // ( 1, 1)
				var h10 = (sbyte) _permutation[(yp + xsb + 2) & 0xFF]; // ( 2, 1)
				yp = _permutation[(ysb + 2) & 0xFF];
				var h11 = (sbyte) _permutation[(yp + xsb + 1) & 0xFF]; // ( 1, 2)
				var h12 = (sbyte) _permutation[(yp + xsb + 2) & 0xFF]; // ( 2, 2)

				value = Kernels(dx, dy, h1, h2, h3,
					h4, h5, h6, h7, h8, h9, h10, h11, h12);
			}
			else
			{
				// We're inside the y>x triangle of the rhombus

				// Get our 12 surrounding vertex values
				short yp = _permutation[(ysb - 1) & 0xFF];
				var h1 = (sbyte) _permutation[(yp + xsb - 1) & 0xFF]; // (-1,-1)
				var h4 = (sbyte) _permutation[(yp + xsb + 0) & 0xFF]; // ( 0,-1)
				yp = _permutation[(ysb + 0) & 0xFF];
				var h2 = (sbyte) _permutation[(yp + xsb - 1) & 0xFF]; // (-1, 0)
				var h5 = (sbyte) _permutation[(yp + xsb + 0) & 0xFF]; // ( 0, 0)
				var h8 = (sbyte) _permutation[(yp + xsb + 1) & 0xFF]; // ( 1, 0)
				yp = _permutation[(ysb + 1) & 0xFF];
				var h3 = (sbyte) _permutation[(yp + xsb - 1) & 0xFF]; // (-1, 1)
				var h6 = (sbyte) _permutation[(yp + xsb + 0) & 0xFF]; // ( 0, 1)
				var h9 = (sbyte) _permutation[(yp + xsb + 1) & 0xFF]; // ( 1, 1)
				var h11 = (sbyte) _permutation[(yp + xsb + 2) & 0xFF]; // ( 2, 1)
				yp = _permutation[(ysb + 2) & 0xFF];
				var h7 = (sbyte) _permutation[(yp + xsb + 0) & 0xFF]; // ( 0, 2)
				var h10 = (sbyte) _permutation[(yp + xsb + 1) & 0xFF]; // ( 1, 2)
				var h12 = (sbyte) _permutation[(yp + xsb + 2) & 0xFF]; // ( 2, 2)

				value = Kernels(dy, dx, h1, h2, h3,
					h4, h5, h6, h7, h8, h9, h10, h11, h12);
			}
			return value/Normalization;
		}

		private static double Kernels(double dx, double dy,
			sbyte h1, sbyte h2, sbyte h3, sbyte h4, sbyte h5, sbyte h6,
			sbyte h7, sbyte h8, sbyte h9, sbyte h10, sbyte h11, sbyte h12)
		{
			double value = 0;

			double dxv = dx + 1 + 2*Stretch;
			double dyv = dy + 1 + 2*Stretch;
			double attn = KernelRadius - dxv*dxv - dyv*dyv;
			if (attn > 0)
			{
				attn *= attn;
				value += attn*attn*h1;
			}

			dxv = dx + 0 + Stretch;
			dyv = dy + 1 + Stretch;
			attn = KernelRadius - dxv*dxv - dyv*dyv;
			// if (attn > 0) {
			attn *= attn;
			value += attn*attn*h2;
			// }

			dxv = dx - 1;
			dyv = dy + 1;
			attn = KernelRadius - dxv*dxv - dyv*dyv;
			if (attn > 0)
			{
				attn *= attn;
				value += attn*attn*h3;
			}

			dxv = dx + 1 + Stretch;
			dyv = dy + 0 + Stretch;
			attn = KernelRadius - dxv*dxv - dyv*dyv;
			if (attn > 0)
			{
				attn *= attn;
				value += attn*attn*h4;
			}

			dxv = dx + 0;
			dyv = dy + 0;
			attn = KernelRadius - dxv*dxv - dyv*dyv;
			// if (attn > 0) {
			attn *= attn;
			value += attn*attn*h5;
			// }

			dxv = dx - 1 - Stretch;
			dyv = dy + 0 - Stretch;
			attn = KernelRadius - dxv*dxv - dyv*dyv;
			// if (attn > 0) {
			attn *= attn;
			value += attn*attn*h6;
			// }

			dxv = dx - 2 - 2*Stretch;
			dyv = dy + 0 - 2*Stretch;
			attn = KernelRadius - dxv*dxv - dyv*dyv;
			if (attn > 0)
			{
				attn *= attn;
				value += attn*attn*h7;
			}

			dxv = dx + 0 - Stretch;
			dyv = dy - 1 - Stretch;
			attn = KernelRadius - dxv*dxv - dyv*dyv;
			// if (attn > 0) {
			attn *= attn;
			value += attn*attn*h8;
			// }

			dxv = dx - 1 - 2*Stretch;
			dyv = dy - 1 - 2*Stretch;
			attn = KernelRadius - dxv*dxv - dyv*dyv;
			// if (attn > 0) {
			attn *= attn;
			value += attn*attn*h9;
			// }

			dxv = dx - 2 - 3*Stretch;
			dyv = dy - 1 - 3*Stretch;
			attn = KernelRadius - dxv*dxv - dyv*dyv;
			// if (attn > 0) {
			attn *= attn;
			value += attn*attn*h10;
			// }

			dxv = dx - 1 - 3*Stretch;
			dyv = dy - 2 - 3*Stretch;
			attn = KernelRadius - dxv*dxv - dyv*dyv;
			if (attn > 0)
			{
				attn *= attn;
				value += attn*attn*h11;
			}

			dxv = dx - 2 - 4*Stretch;
			dyv = dy - 2 - 4*Stretch;
			attn = KernelRadius - dxv*dxv - dyv*dyv;
			if (attn > 0)
			{
				attn *= attn;
				value += attn*attn*h12;
			}

			return value;
		}

		private static int FastFloor(double x)
		{
			var xi = (int) x;
			return x < xi ? xi - 1 : xi;
		}
	}
}