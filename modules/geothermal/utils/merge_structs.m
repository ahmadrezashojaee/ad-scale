function struct = merge_structs(struct1, struct2)
    struct = struct1;
    fn = reshape(fieldnames(struct2), 1, []);
    for n = fn
        struct.(fn{1}) = struct2.(fn{1});
    end
end