using System.Diagnostics.Contracts;

namespace BitFn.NoiseLibrary
{
	public interface INoiseService2D
	{
		[Pure]
		double ComputeNoise(double x, double y);
	}
}