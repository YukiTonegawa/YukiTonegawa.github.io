function activate_copy_button(text_id, button_id) {
    const btn = document.getElementById(button_id);
    const txt = document.getElementById(text_id).value;
    btn.addEventListener('click', () => {
        navigator.clipboard.writeText(txt);
        btn.innerHTML = 'OK';
        setTimeout(() => (btn.innerHTML = 'Copy'), 1000);
    });
}

function cmp_to_string(cmp) {
    let mp = new Map();
    for (let i = 0; i < cmp.length; i++) {
        if (!mp.has(cmp[i])) {
            mp.set(cmp[i], []);
        }
        mp.get(cmp[i]).push(i);
    }
    let result = mp.size + "\n";
    for (const [id, ar] of mp) {
        let cnt = ar.length;
        for (let j = 0; j < cnt; j++) {
            result += (ar[j] + 1) + (j + 1 == cnt ? "\n" : " ");
        }
    }
    return result;
}