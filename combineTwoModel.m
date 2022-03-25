function [mixedModel] = combineTwoModel(model1, model2, addedrxnlist)

%defult addedrxnlist: all diffent reactions but ex reaations
difRxn = setdiff(model2.rxns, model1.rxns);
%difRxn(strmatch('EX_',difRxn))=[]; %remove the EX reactions

if isempty(addedrxnlist)
	addedrxnlist = difRxn;
end	

Addedmodel = model1;
for r = 1:length(addedrxnlist)
	addrxnid = find(strcmp(model2.rxns, addedrxnlist(r)));
	if (addrxnid)
		metsid = find(model2.S(:,addrxnid));
		Srxns = model2.S(metsid,addrxnid);
		addedmet = model2.mets(metsid);
        addrxnlb = model2.lb(addrxnid);
        addrxnub = model2.ub(addrxnid);
		Addedmodel = addReaction(Addedmodel,char(addedrxnlist(r)),'metaboliteList',addedmet,'stoichCoeffList',Srxns, 'lowerBound',addrxnlb,'upperBound', addrxnub);
	else
		addedrxnlist(r)
	end
end

mixedModel = Addedmodel;