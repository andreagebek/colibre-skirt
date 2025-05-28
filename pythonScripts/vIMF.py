import numpy as np

def IMFParameters(snapshot):
    def get_snapshot_param_float(snapshot, param_name: str) -> float:
        try:
            return float(snapshot.metadata.parameters[param_name].decode("utf-8"))
        except KeyError:
            raise KeyError(f"Parameter {param_name} not found in snapshot metadata.")


    try:
        # Extract IMF parameters from snapshot metadata
        alpha_min = get_snapshot_param_float(
            snapshot, "COLIBREFeedback:IMF_HighMass_slope_minimum"
        )
        alpha_max = get_snapshot_param_float(
            snapshot, "COLIBREFeedback:IMF_HighMass_slope_maximum"
        )
        sigma = get_snapshot_param_float(
            snapshot, "COLIBREFeedback:IMF_sigmoid_inverse_width"
        )
        pivot_density = get_snapshot_param_float(
            snapshot, "COLIBREFeedback:IMF_sigmoid_pivot_CGS"
        )

    except:
        alpha_min = alpha_max = -2.3
        sigma = pivot_density = 1

    return (alpha_min, alpha_max, sigma, pivot_density)

def sigmoid(X,x_low,x_high,sigma,pivot):
        return (x_low-x_high)/(1+np.exp(sigma*np.log10(X/pivot)))+x_high;


