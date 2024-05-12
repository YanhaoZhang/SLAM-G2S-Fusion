function vec =  BCH_Residual_rh (eta, x)

    ad_eta = SE3.adHat(eta);
    
    ad_x = SE3.adHat(x);
    
    vec = ad_x * ad_x * eta / 12 ...
             + ad_eta * ad_x * ad_x * eta / 24;
    
end