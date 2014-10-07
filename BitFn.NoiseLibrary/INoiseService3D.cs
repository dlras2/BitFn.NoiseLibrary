using System.Diagnostics.Contracts;

namespace BitFn.NoiseLibrary
{
	public interface INoiseService3D
	{
		[Pure]
		double ComputeNoise(double x, double y, double z);
	}
}