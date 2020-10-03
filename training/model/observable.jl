const observables = [
    "Phosphorylated_Akt"
    "Phosphorylated_ERK"
    "Phosphorylated_cFos"
];

function observables_index(observable_name::String)::Int

    return findfirst(isequal(observable_name),observables)
end