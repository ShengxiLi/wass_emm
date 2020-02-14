% calculate emperial distributions
function pdf_emp = emp_pdf(t_vec, proj_samples)
    [pdf_emp,~] = histcounts(proj_samples, [t_vec t_vec(end)]);
    pdf_emp = pdf_emp./sum(pdf_emp);