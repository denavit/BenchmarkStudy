classdef notional_load
    
    properties
        factorAlwaysAdditive = 0.0;
        factorConditionalB2 = 0.0;
        limitB2 = Inf;
    end
    
    methods
        function obj = notional_load(factorAlwaysAdditive,factorConditionalB2,limitB2)
            if nargin >= 1
                obj.factorAlwaysAdditive = factorAlwaysAdditive;
            end
            if nargin >= 2
                obj.factorConditionalB2 = factorConditionalB2;
            end
            if nargin >=3
                obj.limitB2 = limitB2;
            end
        end
        
        function Hwnl = HwithNotionalLoad(obj,H,Ptotal,B2)
            assert(Ptotal<=0,'Ptotal should be negative, Ptotal = %g',Ptotal);
            if B2 < obj.limitB2
                % Low second-order effects - factorConditionalB2 is minimum
                Hwnl = max(H,obj.factorConditionalB2*-Ptotal) + obj.factorAlwaysAdditive*-Ptotal;
                
            else
                % High second-order effects - factorConditionalB2 is additive
                Hwnl = H + obj.factorConditionalB2*-Ptotal + obj.factorAlwaysAdditive*-Ptotal;
                
            end
        end
        
        function Hest = LoadEstimate(obj,Hest1,Ptotal,B2)
            assert(Ptotal<=0,'Ptotal should be negative, Ptotal = %g',Ptotal);
            Hest = max(Hest1 - obj.factorAlwaysAdditive*-Ptotal,0);
            if  B2 < obj.limitB2
                Hest = max(Hest,1.1*obj.factorConditionalB2*-Ptotal);
            end
            
            % Hest = max(M2/obj.maxFirstOrderMoment(1.0) - notionalLoadRatio*(1+obj.gamma)*(-P),0);
            % if strcmpi(notionalLoadRatioType,'minimum') || (strcmpi(notionalLoadRatioType,'AISC2010 C2.2b(4)') && obj.driftRatio(P) < 1.7)
            %     Hest = max(Hest,1.1*notionalLoadRatio*(1+obj.gamma)*(-P));
            % end
        end
    end
    
end

